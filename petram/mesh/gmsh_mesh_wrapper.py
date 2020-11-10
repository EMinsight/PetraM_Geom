from __future__ import print_function

from petram.debug import timeit, use_profiler
import petram.geom.gmsh_config as gmsh_config
import os
import numpy as np
import gmsh
import time
import tempfile
import traceback
from collections import defaultdict

from scipy.spatial import cKDTree
import multiprocessing as mp
from six.moves.queue import Empty as QueueEmpty

from collections import OrderedDict
Algorithm2D = OrderedDict((("MeshAdap", 1), ("Automatic", 2), ("Delaunay", 5),
                           ("Frontal", 6), ("BAMG", 7), 
                           ("FrrontalQuad", 8), ("Paking of parallelograms", 9),
                           ("default", 2)))
Algorithm3D = OrderedDict((("Delaunay", 1), 
                           ("Frontal", 4),
                           ("HXT", 10), ("MMG3D", 7),
                           ("R-tree", 9), ("default", 1)))
AlgorithmR = OrderedDict((("Simple", 0), ("Blossom", 1),
                          ("SimpleFullQuad", 2), ("BlossomFullQuad", 3),
                          ("default", 3),))
HighOrderOptimize = OrderedDict((("none", 0),
                                 ("optimization", 1),
                                 ("elastic+optimization", 2),
                                 ("elastic", 3),
                                 ("fast curving", 4)))
debug = True
debug2 = False

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GMSHMeshWrapper')

def dprint(*args):
    if debug:
        print(*args)


def dprint2(*args):
    if debug2:
        print(*args)


def get_vertex_geom_zie():
    from collections import defaultdict

    lcar = defaultdict(lambda: np.inf)

    for dim, tag in gmsh.model.getEntities(1):
        x1, y1, z1, x2, y2, z2 = gmsh.model.getBoundingBox(dim, tag)
        s = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
        bdimtags = gmsh.model.getBoundary(((dim, tag,),), oriented=False)
        for bdim, btag in bdimtags:
            lcar[btag] = min((lcar[btag], s))
    return dict(lcar)


def process_text_tags(dim=1, check=True):
    '''
    convert text tags input to dimtags
    tags = '1,2,3', 'all', or 'remaining'
    '''
    def func2(method):
        def method2(self, done, params, tags, *args, **kwargs):
            if tags == 'remaining':
                dimtags = gmsh.model.getEntities(dim)
                if check:
                    dimtags = [(dim, x)
                               for xx, x in dimtags if not x in done[dim]]
                else:
                    pass
            elif tags == 'all':
                dimtags = gmsh.model.getEntities(dim)
            else:
                tags = [int(x) for x in tags.split(',')]
                if check:
                    dimtags = [(dim, x) for x in tags if not x in done[dim]]
                else:
                    dimtags = [(dim, x) for x in tags]

            return method(self, done, params, dimtags, *args, **kwargs)
        return method2
    return func2


def process_text_tags_sd(dim=1):
    '''
    convert two text tags input to dimtags
    tags = '1,2,3'
    '''
    def func2(method):
        def method2(self, done, params, tags, tags2, *args, **kwargs):
            if tags == 'remaining':
                assert False, "S-D mapping does not support remaining"
            elif tags == 'all':
                assert False, "S-D mapping does not support all"
            else:
                tags = [int(x) for x in tags.split(',')]
                dimtags = [(dim, x) for x in tags]

            if tags2 == 'remaining':
                assert False, "S-D mapping does not support remaining"
            elif tags2 == 'all':
                assert False, "S-D mapping does not support all"
            else:
                tags2 = [int(x) for x in tags2.split(',')]
                dimtags2 = [(dim, x) for x in tags2]
            return method(self, done, params, dimtags, dimtags2, *args, **kwargs)
        return method2
    return func2


def process_text_tags_vsd(dim=3):
    '''
    convert three text tags input to dimtags
    tags = '1,2,3'
    '''
    def func2(method):
        def method2(self, done, params, vtags, tags, tags2, *args, **kwargs):
            vtags = [int(x) for x in vtags.split(',')]
            vdimtags = [(dim, x) for x in vtags]

            tags = [int(x) for x in tags.split(',')]
            dimtags = [(dim-1, x) for x in tags]

            tags2 = [int(x) for x in tags2.split(',')]
            dimtags2 = [(dim-1, x) for x in tags2]
            return method(self, done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)
        return method2
    return func2


def set_restore_maxmin_cl(method):
    def wrapped2(self, *args, **kwargs):
        maxsize = kwargs.get("maxsize", self.clmax)
        minsize = kwargs.get("minsize", self.clmin)

        if maxsize != 0:
            value = maxsize
        else:
            value = self.clmax
        kwargs['maxsize'] = value
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", value)
        
        if minsize != 0:
            value = minsize
        else:
            value = self.clmin
        kwargs['minsize'] = value
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", value)        

        ret = method(self, *args, **kwargs)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.clmax)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.clmin)
        return ret
    return wrapped2

def check_line_orientation(ltag, vtags, pcoord):
    p1 = np.array(gmsh.model.getValue(0, vtags[0], [0]))
    p2 = np.array(gmsh.model.getValue(0, vtags[1], [0]))
    p3 = np.array(gmsh.model.getValue(1, ltag, [pcoord]))

    d1 = np.sqrt(np.sum((p1-p3)**2))
    d2 = np.sqrt(np.sum((p2-p3)**2))

    if d1 > 1e-10 and d2 > 1e-10:
        print("Line endes does not agree with Vertex")
        return 0
    if d1 < d2:
        return 1
    if d2 < d1:
        return 2


def get_nodes_elements(ents, normalize=False):
    mdata = []
    for dim, tag in ents:

        if normalize:
            ndata = gmsh.model.mesh.getNodes(dim, tag, includeBoundary=True)
            pc = np.array(ndata[2])
            ndata = ndata[0], ndata[1], (pc[:-2]-pc[-2])/(pc[-1]-pc[-2])
        else:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
        edata = gmsh.model.mesh.getElements(dim, tag)
        mdata.append((dim, tag, ndata, edata))
    return mdata


class GMSHMeshWrapper():
    workspace_base = 1
    gmsh_init = False

    def __init__(self, meshformat=2.2,
                 CharacteristicLengthMax=1e20,
                 CharacteristicLengthMin=0,
                 EdgeResolution=3,
                 MeshAlgorithm="Automatic",
                 MeshAlgorithm3D="Delaunay",
                 MaxThreads=[1, 1, 1, 1],
                 **kwargs):

        self.queue = kwargs.pop("queue", None)
        self.use_profiler = kwargs.pop("use_profiler", True)
        self.use_expert_mode = kwargs.pop("use_expert_mode", False)
        self.gen_all_phys_entity = kwargs.pop("gen_all_phys_entity", False)
        self.trash = kwargs.pop("trash", '')
        self.edge_tss = kwargs.pop("edge_tss", None)
        self.mesh_sequence = kwargs.pop("mesh_sequence", [])
        self.use_ho = kwargs.pop("use_ho", False)
        self.ho_order = kwargs.pop("ho_order", 2)        
        self.optimize_ho = kwargs.pop("optimize_ho", 0)
        self.optimize_dom = kwargs.pop("optimize_dom", "all")        
        self.mapper_tol = kwargs.pop("mapper_tol", 1e-5)
        
        gmsh.clear()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", meshformat)
        gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)
        gmsh.option.setNumber("Mesh.IgnorePeriodicity", 1)

        gmsh_init = True
        self.add_model('main1')

        # default clmax
        self.clmax = CharacteristicLengthMax
        self.clmin = CharacteristicLengthMin

        self.res = EdgeResolution
        self.algorithm = MeshAlgorithm
        self.algorithm3d = MeshAlgorithm3D
        self.algorithmr = kwargs.pop("MeshAlgorithmR", "default")
        self.maxthreads = MaxThreads    # general, 1D, 2D, 3D: defualt = 1,1,1,1
        self._new_brep = True
        self._name = "GMSH_Mesher"

    def name(self):
        return self._name

    @use_profiler
    def generate(self, brep_input, msh_file, dim=3, finalize=False):
        '''
        generate mesh based on  meshing job sequence.
        brep must be loaed 
        '''
        print("brep input", brep_input)
        if brep_input != '':
            self.load_brep(brep_input)
        if not self._new_brep:
            assert False, "new BREP must be loaded"

        self._new_brep = False
        GMSHMeshWrapper.workspace_base = 1
        self.switch_model('main1')

        gmsh.option.setNumber("Mesh.Algorithm",
                              Algorithm2D[self.algorithm])
        gmsh.option.setNumber("Mesh.Algorithm3D",
                              Algorithm3D[self.algorithm3d])
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm",
                              AlgorithmR[self.algorithmr])
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.clmax)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.clmin)
        
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        gmsh.option.setNumber("General.ExpertMode",
                              1 if self.use_expert_mode else 0)
        self.maxthreads = (3,3,3,3)
        gmsh.option.setNumber("General.NumThreads",   self.maxthreads[0])
        gmsh.option.setNumber("Mesh.MaxNumThreads1D", self.maxthreads[1])
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", self.maxthreads[2])
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", self.maxthreads[3])
        gmsh.option.setNumber("Mesh.Optimize", 0)
        #gmsh.option.setNumber('Geometry.ReparamOnFaceRobust', 1)

        #if self.use_ho:
        #    gmsh.option.setNumber("Mesh.ElementOrder", self.ho_order)
        #    gmsh.option.setNumber("Mesh.HighOrderOptimize",
        #                          HighOrderOptimize[self.optimize_ho])
        
        self.target_entities0 = (gmsh.model.getEntities(3),
                                 gmsh.model.getEntities(2),
                                 gmsh.model.getEntities(1),
                                 gmsh.model.getEntities(0))
        self.target_entities = gmsh.model.getEntities()
        #
        self.vertex_geom_size = get_vertex_geom_zie()

        # set default vertex mesh size
        self.cl_data = {}
        sizes = []
        for tag in self.vertex_geom_size.keys():
            size = self.vertex_geom_size[tag]/self.res
            if size > self.clmax:
                size = self.clmax
            if size <= self.clmin:
                size = self.clmin
                
            self.setCL(((0, tag),), size)
            sizes.append(size)
            #print("Default Point Size", (0, tag), size)
        if len(sizes) > 0:
            print("Default Point Size Max/Min", max(sizes), min(sizes))

        done = [[], [], [], []]
        params = [None]*len(self.mesh_sequence)

        maxdim = max([x for x, tag in gmsh.model.getEntities()])
        maxdim = min([maxdim, dim])
        #adim, idx_dim = self.check_algorith_dim()
        #adim  = self.check_algorith_dim()

        for mdim in range(maxdim+1):
            for idx, sq in enumerate(self.mesh_sequence):
                '''
                if (mdim == adim and
                    idx_dim[mdim] == idx and
                    self.use_ho):
                    #print("turn on high order", self.ho_order)
                    gmsh.option.setNumber("Mesh.ElementOrder", self.ho_order)
                    gmsh.option.setNumber("Mesh.HighOrderOptimize",
                                          HighOrderOptimize[self.optimize_ho])
                '''
                proc, args, kwargs = sq
                f = getattr(self, proc+"_"+str(mdim)+"D")

                if self.queue is not None:
                    self.queue.put((False,
                                    "Processing " + proc+"_"+str(mdim)+"D"))

                if mdim == 3 and idx == len(self.mesh_sequence)-1:
                    gmsh.option.setNumber("Mesh.Optimize", 1)

                done, params[idx] = f(done, params[idx],
                                      *args, **kwargs)
                # gmsh.model.mesh.removeDuplicateNodes()

            for i in range(mdim+1, 4):
                done[i] = []
                
        if self.use_ho:
            # using this option makes "computing connectivity and bad
            # elements very slow"
            gmsh.option.setNumber("Mesh.HighOrderDistCAD", 0)            
            gmsh.option.setNumber("Mesh.ElementOrder", self.ho_order)
            #gmsh.option.setNumber("Mesh.HighOrderOptimize",
            #                       HighOrderOptimize[self.optimize_ho])
            self.hide_all()
            gmsh.model.mesh.generate(3)

            #
            #self.show_all()
            maxdim = max([x for x, tag in gmsh.model.getEntities()])
            if self.optimize_dom.lower() == 'all':
                dimTags = gmsh.model.getEntities(maxdim)
            else:
                dimTags = [(maxdim, int(x)) for x in self.optimize_dom.split(',')]

            #gmsh.option.setNumber("Mesh.MeshOnlyVisible", 0)
            #self.show_all()
            # it is not well-written but I have to include all 
            # boundaries
            # elastic seems to apply for everything anyway ???
            self.show_only(dimTags, recursive=True)            
            if maxdim == 3:
                dimTags1 = gmsh.model.getBoundary(dimTags)
                dimTags2 = gmsh.model.getBoundary(dimTags1)
                dimTags3 = gmsh.model.getBoundary(dimTags2)
                dimTags = list(set(dimTags + dimTags1 + dimTags2 + dimTags3))
                do_ho = True
            elif maxdim == 2:
                dimTags1 = gmsh.model.getBoundary(dimTags)
                dimTags2 = gmsh.model.getBoundary(dimTags1)
                dimTags = list(set(dimTags + dimTags1 + dimTags2))
                do_ho = True
            else:
                do_ho = False
            #do_ho = (do_ho and HighOrderOptimize[self.optimize_ho] != 0)

            if do_ho:
                '''
                if maxdim == 3:
                    gmsh.model.mesh.optimize("Relocate3D", dimTags=[])
                elif maxdim == 2:
                    gmsh.model.mesh.optimize("Relocate2D", dimTags=[])
                else:
                    pass
                '''
                gmsh.option.setNumber("Mesh.HighOrderThresholdMax", 2)
                gmsh.option.setNumber("Mesh.HighOrderThresholdMin", 0.1)

                if HighOrderOptimize[self.optimize_ho] == 1:
                     gmsh.model.mesh.optimize("HighOrder", dimTags=dimTags)                
                elif HighOrderOptimize[self.optimize_ho] == 2:
                     gmsh.model.mesh.optimize("HighOrderElastic", dimTags=dimTags)                
                     gmsh.model.mesh.optimize("HighOrder", dimTags=dimTags)                
                elif HighOrderOptimize[self.optimize_ho] == 3:
                     gmsh.model.mesh.optimize("HighOrderElastic", dimTags=dimTags)
                elif HighOrderOptimize[self.optimize_ho] == 4:
                     gmsh.model.mesh.optimize("HighOrderFastCurving", dimTags=dimTags)
                else:
                    pass

            
        gmsh.model.mesh.removeDuplicateNodes()

        # somehow add_physical is very slow when there are too many physicals...
        # we process the mesh file w/o physical

        use_add_physical = False
        path = os.path.dirname(msh_file)

        if finalize:
            self.queue.put((False,
                            "Adding Physical Entities "))

            if use_add_physical:
                self.add_sequential_physicals()
                gmsh.write(msh_file)
            else:
                print("creating temporary mesh file")
                tmp0 = os.path.join(self.trash, 'tmp0.msh')
                gmsh.write(tmp0)

                print("generating final mesh file", msh_file)
                self.edit_msh_to_add_sequential_physicals(tmp0, msh_file)
                self.edit_msh_to_add_sequential_physicals(tmp0, msh_file+'2.msh')
        else:
            print("creating temporary mesh file")
            tmp0 = os.path.join(self.trash, 'tmp0.msh')
            gmsh.write(tmp0)
            msh_file = tmp0

        return maxdim, done, msh_file

    def load_brep(self, filename):
        self.switch_model('main1')
        self.target = filename

        gmsh.model.occ.importShapes(filename, highestDimOnly=False)
        gmsh.model.occ.synchronize()

        self.geom_info = self.read_geom_info()
        self._new_brep = True
        return self.current

    def read_geom_info(self, dimtags=None):
        from petram.geom.read_gmsh import read_loops3

        self.hide_all()
        gmsh.model.mesh.generate(1)
        # gmsh.model.mesh.removeDuplicateNodes()

        #ptx, p, l, s, v = read_loops2(gmsh)
        return read_loops3(gmsh, dimtags)

    def add_model(self, name):
        gmsh.model.add(name)
        name = gmsh.model.list()[-1]
        self.current = name
        return name

    def switch_model(self, name='main1'):
        gmsh.model.setCurrent(name)
        self.current = name
        return name

    def prep_workspace(self):
        name = 'ws'+str(GMSHMeshWrapper.workspace_base)
        GMSHMeshWrapper.workspace_base += 1

        self.add_model(name)
        gmsh.model.setCurrent(name)
        gmsh.model.occ.importShapes(self.target, highestDimOnly=False)
        self.hide_all()

        gmsh.model.mesh.generate(1)
        gmsh.model.occ.synchronize()

        return name

    def show_only(self, dimtags, recursive=False):
        self.hide_all()
        gmsh.model.setVisibility(dimtags, True, recursive=recursive)

        ent = gmsh.model.getEntities()

        # we hide the surfaces genrated from virtual operation always.
        if self.current == 'main1':
            vent = list(set(ent).difference(self.target_entities))
            #vent = [x for x in ent if not x in self.target_entities]
            #print("hiding virtual", vent)
            if len(vent) > 0:
                gmsh.model.setVisibility(vent, False, recursive=True)

    def hide(self, dimtags, recursive=False):
        gmsh.model.setVisibility(dimtags, False, recursive=recursive)

    def hide_all(self):
        ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, False)

    def show_all(self):
        if self.current == 'main1':
            ent = self.target_entities
        else:
            ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, True)

    @timeit
    def delete_all_except_face(self, tags):
        dimtags = [(2, x) for x in tags]
        ret_3D = gmsh.model.getEntities(3)

        del_list = []
        rem3D_list = []
        rem2D_list = []
        for x in ret_3D:
            bdrs = gmsh.model.getBoundary(x, combined=False, oriented=False)
            print(x, bdrs)
            print(set(bdrs).intersection(dimtags))
            if len(set(bdrs).intersection(dimtags)) == 0:
                print('delete ', x)
                del_list.append(x)
            else:
                rem3D_list.append(x)
                rem2D_list.extend(bdrs)
        print("delete this (3D)#", len(del_list))
        gmsh.model.occ.remove(del_list, recursive=True)

        print("delete bodies #", len(rem3D_list))
        gmsh.model.occ.remove(rem3D_list, recursive=False)

        del_list = list(set(rem2D_list).difference(dimtags))
        print("delete faces #", len(del_list))
        gmsh.model.occ.remove(del_list, recursive=True)

        print('calling this')
        gmsh.model.occ.synchronize()
        print('done')
    '''
    def transfer_mesh(self, dimtags, ws1, ws2, resursive = True):
        if resursive:
            dimtags = self.expand_dimtags(dimtags)

        old_model = self.current
        
        data = []
        self.switch_model(ws1)
        for dim, tag in dimtags:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            data.append((dim, tag, ndata, edata))
        self.switch_model(ws2)            
        for dim, tag, ndata, edata in data:
            gmsh.model.mesh.setNodes(dim, tag, ndata[0], ndata[1], ndata[2])
            gmsh.model.mesh.setElements(dim, tag, edata[0], edata[1], edata[2])

        self.switch_model(old_model)
    '''

    def save_mesh(self, filename, ws=''):
        '''
        save mesh
        '''
        current = self.current
        self.switch_model(ws)
        gmsh.write(filename)
        self.switch_model(current)

    def save_mesh_all(self, filename_base):
        '''
        save mesh for all models
        '''
        workspaces = gmsh.model.list()
        for ws in workspaces:
            f = filename_base+'_'+ws+'.msh'
            self.save_mesh(f, ws=ws)

    def expand_dimtags(self, dimtags, return_dim):
        '''
        expand all subdims
        '''
        while dimtags[0][0] != return_dim:
            bdimtags = []
            for dimtag in dimtags:
                bdimtags.extend(gmsh.model.getBoundary([dimtag,],
                                                      combined=False,
                                                      oriented=False))
            dimtags = tuple(set(tuple(bdimtags)))
        return list(dimtags)
        '''
        volumes = []
        surfaces = []
        lines = []
        points = []
        
        for dim, tag in dimtags:
            if dim == 3:
                volumes.append(tag)
            if dim == 2:
                surfaces.append(tag)
            if dim == 1:
                lines.append(tag)
            if dim == 0:
                points.append(tag)

        for v in volumes:
            surfaces.extend(self.geom_info[4][v])

        surfaces = list(set(surfaces))
        for s in surfaces:
            lines.extend(self.geom_info[3][s])

        lines = list(set(lines))
        for l in lines:
            points.extend(self.geom_info[2][l])

        points = list(set(points))

        ans = ([(0, p) for p in points] +
               [(1, l) for l in lines] +
               [(2, s) for s in surfaces] +
               [(3, v) for v in volumes])

        if return_dim is not None:
            return [(return_dim, x[1]) for x in ans if x[0] == return_dim]
        return ans
        '''

    def list_entities(self):
        for x in gmsh.model.list():
            gmsh.model.setCurrent(x)
            #print(x, gmsh.model.getEntities())
        gmsh.model.setCurrent(self.current)

    def edit_msh_to_add_sequential_physicals(self, tmp_file, filename, verbose=True):

        from petram.geom.read_gmsh import gmsh_element_type, gmsh_element_dim

        lines = OrderedDict()

        fid = open(tmp_file, 'r')

        def readline1():
            ret = fid.readline()
            new_lines1.append(ret)

        l = fid.readline()
        while l:
            line = l.strip()
            if len(line) == 0:
                l = fid.readline()
                continue

            if line[:4] == '$End':
                pass
            elif line[0] == '$':
                section = line[1:]
                lines[section] = []
            else:
                lines[section].append(line)
            l = fid.readline()
        fid.close()
        
        # process element
        ndims = [defaultdict(int),
                 defaultdict(int),
                 defaultdict(int),
                 defaultdict(int)]

        elines = [[], [], [], []]

        for k, l in enumerate(lines["Elements"][1:]):
            xx = l.split(' ')
            el_type = int(xx[1])
            el_num = int(xx[4])
            xx[3] = str(el_num)
            dd = gmsh_element_dim[el_type]
            ndims[dd][el_num] = ndims[dd][el_num] + 1

            elines[dd].append(' '.join(xx))

        if not self.gen_all_phys_entity and len(ndims[3]) != 0:
            elines2 = elines[2] + elines[3]
            nphys = len(ndims[2]) + len(ndims[3])
            ndims[0] = {}
            ndims[1] = {}
            if verbose:
                print("Adding " + str(len(ndims[2])) + " Surface(s)")
                print("Adding " + str(len(ndims[3])) + " Volume(s)")

        elif not self.gen_all_phys_entity and len(ndims[2]) != 0:
            elines2 = elines[1] + elines[2]
            nphys = len(ndims[1]) + len(ndims[2])
            ndims[0] = {}
            if verbose:
                print("Adding " + str(len(ndims[1])) + " Line(s)")
                print("Adding " + str(len(ndims[2])) + " Surface(s)")

        else:
            elines2 = elines[0] + elines[1] + elines[2] + elines[3]
            nphys = len(ndims[0]) + len(ndims[1]) + \
                len(ndims[2]) + len(ndims[3])
            if verbose:
                print("Adding " + str(len(ndims[0])) + " Point(s)")
                print("Adding " + str(len(ndims[1])) + " Line(s)")
                print("Adding " + str(len(ndims[2])) + " Surfac(s)")
                print("Adding " + str(len(ndims[3])) + " Volume(s)")

        # renumber elements
        elines3 = []
        for k, l in enumerate(elines2):
            xx = l.split(' ')
            elines3.append(' '.join([str(k+1)] + xx[1:]))
        elines3 = [str(len(elines3))]+elines3
        lines["Elements"] = elines3

        phys_names = [str(nphys)]
        for l in list(ndims[0]):
            phys_names.append(" ".join(["0", str(l), '"point'+str(l)+'"']))
        for l in list(ndims[1]):
            phys_names.append(" ".join(["1", str(l), '"line'+str(l)+'"']))
        for l in list(ndims[2]):
            phys_names.append(" ".join(["2", str(l), '"surface'+str(l)+'"']))
        for l in list(ndims[3]):
            phys_names.append(" ".join(["3", str(l), '"volume'+str(l)+'"']))

        lines["PhysicalNames"] = phys_names

        def write_section(fid, lines, sec):
            fid.write("$"+sec+"\n")
            fid.write("\n".join(lines[sec])+"\n")
            fid.write("$End"+sec+"\n")

        rest_sec = [x for x in list(lines) if x !=
                    "MeshFormat" and x != "PhysicalNames"]

        fid = open(filename, "w")

        write_section(fid, lines, "MeshFormat")
        write_section(fid, lines, "PhysicalNames")
        for sec in rest_sec:
            write_section(fid, lines, sec)

        fid.close()
        #from shutil import copyfile
        #copyfile(filename, filename+'.bk')

    def add_sequential_physicals(self, verbose=True):
        '''
        add sequencial physical entity numbers
        '''
        if self.queue is not None:
            self.queue.put((False, "Adding Physicals..."))

        ent = gmsh.model.getEntities(dim=3)
        max_dim = 0

        if len(ent) > 0:
            max_dim = 3
            if verbose:
                print("Adding " + str(len(ent)) + " Volume(s)")
        for k, x in enumerate(ent):
            # if len(gmsh.model.getPhysicalGroupsForEntity(3, x[1])) > 0:
            #    continue
            value = gmsh.model.addPhysicalGroup(3, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(3, value, 'volume'+str(value))

        ent = gmsh.model.getEntities(dim=2)
        if len(ent) > 0:
            if max_dim == 0:
                max_dim = 2
            if verbose:
                print("Adding " + str(len(ent)) + " Surface(s)")

        for k, x in enumerate(ent):
            # if len(gmsh.model.getPhysicalGroupsForEntity(2, x[1])) > 0:
            #    continue
            value = gmsh.model.addPhysicalGroup(2, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(2, value, 'surface'+str(value))

        if not self.gen_all_phys_entity and max_dim == 3:
            return
        ent = gmsh.model.getEntities(dim=1)
        if len(ent) > 0:
            if max_dim == 0:
                max_dim = 1
            if verbose:
                print("Adding " + str(len(ent)) + " Line(s)")

        for k, x in enumerate(ent):
            # if len(gmsh.model.getPhysicalGroupsForEntity(1, x[1])) > 0:
            #    continue
            value = gmsh.model.addPhysicalGroup(1, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(1, value, 'line'+str(value))

        if not self.gen_all_phys_entity and max_dim == 2:
            return
        ent = gmsh.model.getEntities(dim=0)
        if len(ent) > 0:
            if verbose:
                print("Adding " + str(len(ent)) + " Point(s)")
        for k, x in enumerate(ent):
            # if len(gmsh.model.getPhysicalGroupsForEntity(0, x[1])) > 0:
            #    continue
            value = gmsh.model.addPhysicalGroup(0, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(0, value, 'point'+str(value))

    def merge_text(self, geo_text):
        handle, geo_filename = tempfile.mkstemp(suffix='.geo')
        text = geo_text.encode()
        #print("writing this", text)
        os.write(handle, text)
        os.close(handle)
        gmsh.merge(geo_filename)
        os.remove(geo_filename)

    def check_algorith_dim(self):
        dims = {'cl': 0,
                'freevolume':3,
                'freeface':2,
                'freeedge':1,
                'transfinite_volume':3,
                'transfinite_face':2,
                'transfinite_edge':1,
                'recombine_surface':0,
                'copyface':2,
                'extrude_face':2,
                'revolve_face':2,
                'mergetxt':0,
                }
        d = []
        for sq in self.mesh_sequence:
            proc, _args, _kwargs = sq
            d.append(dims[proc])

        return max(d)

    def setCL(self, dimtags, size):
        for dimtag in dimtags:
           self.cl_data[dimtag] = size
        gmsh.model.mesh.setSize(dimtags, size)        
    '''
    Low-level implementation at each mesh dim
    '''
    @process_text_tags(dim=0)
    def cl_0D(self, done, params, dimtags, size, overwrite=True):
        '''
        Set Charcteristic length of points
        By default, it overwrite whatever it is set
        '''

        if not overwrite:
            dimtags = [(0, x) for x in tags]
        else:
            dimtags = [x for x in dimtags if not x[0] in done[0]]
        if len(dimtags) == 0:
            return
        self.show_only(dimtags)
        self.setCL(dimtags, size)

        for dim, tag in dimtags:
            if not tag in done[0]:
                done[0].append(tag)
        gmsh.model.mesh.generate(0)
        return done, params

    @process_text_tags(dim=0)
    def cl_1D(self, done, params, *args, **kwargs):
        return done, params

    @process_text_tags(dim=0)
    def cl_2D(self, done, params, *args, **kwargs):
        return done, params

    @process_text_tags(dim=0)
    def cl_3D(self, done, params, *args, **kwargs):
        return done, params

    @process_text_tags(dim=2, check=False)
    def recombine_surface_0D(self, done, params, dimtags):
        for dim, tag in dimtags:
            if dim != 2:
                continue
            gmsh.model.mesh.setRecombine(2, tag)
        return done, params

    @process_text_tags(dim=2, check=False)
    def recombine_surface_1D(self, done, params, dimtags):
        return done, params

    @process_text_tags(dim=2, check=False)
    def recombine_surface_2D(self, done, params, dimtags):
        return done, params

    @process_text_tags(dim=2, check=False)
    def recombine_surface_3D(self, done, params, dimtags):
        return done, params

    # freevolume
    @set_restore_maxmin_cl
    @process_text_tags(dim=3)
    def freevolume_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize = kwargs.pop("maxsize", 1e20)
        minsize = kwargs.pop("minsize", 0.0)
        res = kwargs.pop("resolution",  np.inf)

        done[3].extend([x for dim, x in dimtags])

        embeds = [x for x in kwargs.pop("embed_s",  '').split(',')]
        embeds = [int(x) for x in embeds if len(x) > 0]
        if len(embeds) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(2, embeds, dimtags[0][0], dimtags[0][1])

        embedl = [x for x in kwargs.pop("embed_l",  '').split(',')]
        embedl = [int(x) for x in embedl if len(x) > 0]
        if len(embedl) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(1, embedl, dimtags[0][0], dimtags[0][1])

        embedp = [x for x in kwargs.pop("embed_p",  '').split(',')]
        embedp = [int(x) for x in embedp if len(x) > 0]
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(0, embedp, dimtags[0][0], dimtags[0][1])

        dimtags.extend([(2, x) for x in embeds])
        sdimtags = self.expand_dimtags(dimtags, return_dim=2)
        done[2].extend([x for dim, x in sdimtags if not x in done[2]])

        dimtags.extend([(1, x) for x in embedl])
        ldimtags = self.expand_dimtags(dimtags, return_dim=1)
        done[1].extend([x for dim, x in ldimtags if not x in done[1]])

        dimtags.extend([(0, x) for x in embedp])
        dimtags = self.expand_dimtags(dimtags, return_dim=0)

        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize:
                size = maxsize
            if size < minsize:
                size = minsize
            self.setCL(((0, tag), ), size)                
            #print("Volume Set Point Size", (0, tag), size)
            done[0].append(tag)
        gmsh.model.mesh.generate(0)
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=3)
    def freevolume_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[3].extend([x for dim, x in dimtags])

        embeds = [x for x in kwargs.pop("embed_s",  '').split(',')]
        embeds = [(2, int(x)) for x in embeds if len(x) > 0]
        embedl = [x for x in kwargs.pop("embed_l",  '').split(',')]
        embedl = [(1, int(x)) for x in embedl if len(x) > 0]

        dimtags.extend(embeds)
        sdimtags = self.expand_dimtags(dimtags, return_dim=2)
        done[2].extend([x for dim, x in sdimtags if not x in done[2]])

        dimtags.extend(embedl)

        dimtags = self.expand_dimtags(dimtags, return_dim=1)

        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]
        #tags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]

        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in dimtags])
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=3)
    def freevolume_2D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[3].extend([x for dim, x in dimtags])

        embeds = [x for x in kwargs.pop("embed_s",  '').split(',')]
        embeds = [(2, int(x)) for x in embeds if len(x) > 0]
        dimtags.extend(embeds)

        dimtags = self.expand_dimtags(dimtags, return_dim=2)

        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[2]]
        #tags = [(dim, tag) for dim, tag in dimtags if not tag in done[2]]

        self.show_only(dimtags)
        #print("2D meshing for ", dimtags)
        alg2d = kwargs.get("alg2d", "default")
        if alg2d != 'default':
            gmsh.option.setNumber("Mesh.Algorithm",
                                  Algorithm2D[alg2d])
            gmsh.model.mesh.generate(2)
            gmsh.option.setNumber("Mesh.Algorithm",
                                  Algorithm2D[self.algorithm])
        else:
            gmsh.model.mesh.generate(2)        
        done[2].extend([x for dim, x in dimtags])
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=3)
    def freevolume_3D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags if not tag in done[3]]
        self.show_only(tags, recursive=True)

        alg3d = kwargs.get("alg3d", "default")
        if alg3d != 'default':
            gmsh.option.setNumber("Mesh.Algorithm3D",
                                  Algorithm3D[alg3d])
            gmsh.model.mesh.generate(3)
            gmsh.option.setNumber("Mesh.Algorithm3D",
                                  Algorithm3D[self.algorithm3d])
        else:
            gmsh.model.mesh.generate(3)            
        
        done[3].extend([x for dim, x in tags])
        return done, params

    # freeface
    @set_restore_maxmin_cl
    @process_text_tags(dim=2)
    def freeface_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize = kwargs.pop("maxsize", 1e20)
        minsize = kwargs.pop("minsize", 0.0)
        res = kwargs.pop("resolution", np.inf)

        done[2].extend([x for dim, x in dimtags])

        embedl = [x for x in kwargs.pop("embed_l",  '').split(',')]
        embedl = [int(x) for x in embedl if len(x) > 0]
        if len(embedl) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(1, embedl, dimtags[0][0], dimtags[0][1])

        embedp = [x for x in kwargs.pop("embed_p",  '').split(',')]
        embedp = [int(x) for x in embedp if len(x) > 0]
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(0, embedp, dimtags[0][0], dimtags[0][1])

        dimtags.extend([(1, x) for x in embedl])
        dimtags.extend([(0, x) for x in embedp])
        dimtags = self.expand_dimtags(dimtags, return_dim=0)

        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize:
                size = maxsize
            if size < minsize:
                size = minsize
            self.setCL(((0, tag), ), size)
            #print("Face Set Point Size", (0, tag), size)
            done[0].append(tag)
        gmsh.model.mesh.generate(0)
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=2)
    def freeface_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[2].extend([x for dim, x in dimtags])

        dimtags = self.expand_dimtags(dimtags, return_dim=1)
        embedl = [x for x in kwargs.pop("embed_l",  '').split(',')]
        embedl = [(1, int(x)) for x in embedl if len(x) > 0]
        dimtags.extend(embedl)

        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]
        #tags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]
        self.show_only(dimtags)

        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 1)
        gmsh.model.mesh.generate(1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
        done[1].extend([x for dim, x in dimtags])

        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=2)
    def freeface_2D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags]
        self.show_only(dimtags)

        alg2d = kwargs.get("alg2d", "default")
        if alg2d != 'default':
            gmsh.option.setNumber("Mesh.Algorithm",
                                    Algorithm2D[alg2d])
            gmsh.model.mesh.generate(2)
            gmsh.option.setNumber("Mesh.Algorithm",
                                   Algorithm2D[self.algorithm])
        else:
            gmsh.model.mesh.generate(2)            
        
        done[2].extend([x for dim, x in tags])
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=2)
    def freeface_3D(self, done, params, tags, *args, **kwargs):
        return done, params

    # freeedge
    @set_restore_maxmin_cl
    @process_text_tags(dim=1)
    def freeedge_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize = kwargs.pop("maxsize", 1e20)
        minsize = kwargs.pop("minsize", 0.0)
        res = kwargs.pop("resolution", np.inf)

        embedp = [x for x in kwargs.pop("embed_p",  '').split(',')]
        embedp = [int(x) for x in embedp if len(x) > 0]
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            gmsh.model.mesh.embed(0, embedp, dimtags[0][0], dimtags[0][1])

        done[1].extend([x for dim, x in dimtags])

        dimtags = self.expand_dimtags(dimtags, return_dim=0)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize:
                size = maxsize
            if size <= minsize:
                size = minsize
            self.setCL(((0, tag), ), size)                
            done[0].append(tag)
        gmsh.model.mesh.generate(0)
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=1)
    def freeedge_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in tags])
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=1)
    def freeedge_2D(self, done, params, tags, *args, **kwargs):
        return done, params

    @set_restore_maxmin_cl
    @process_text_tags(dim=1)
    def freeedge_3D(self,  done, params, tags, *args, **kwargs):
        return done, params

    # transfinite_edge
    @process_text_tags(dim=1)
    def transfinite_edge_0D(self, done, params, dimtags, *args, **kwargs):
        #meher.add('transfinite_line', gid, nseg=nseg, progression = p,  bump = b)
        nseg = kwargs.get('nseg', 100) + 1

        meshType = 'Progression'
        coef = kwargs.get('progression', 1.0)

        bump = kwargs.get('bump', 1)
        if bump != 1:
            meshType = 'Bump'
            coef = bump

        for dim, tag in dimtags:
            if tag in done[1]:
                continue
            gmsh.model.mesh.setTransfiniteCurve(tag, nseg,
                                                meshType=meshType,
                                                coef=coef)
            done[1].append(tag)
        return done, params

    @process_text_tags(dim=1)
    def transfinite_edge_1D(self, done, params, dimtags, *args, **kwargs):
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in dimtags])
        return done, params

    def transfinite_edge_2D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    def transfinite_edge_3D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    # transfinite_face
    @process_text_tags(dim=2)
    def transfinite_surface_0D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    @process_text_tags(dim=2)
    def transfinite_surface_1D(self, done, params, dimtags, *args, **kwargs):
        arrangement = kwargs.get('arrangement', 'Left')
        cornerTags = kwargs.get('corner', [])
        #print('Corner', cornerTags)

        # for now, we don't do anything
        # we could add a step to try trasnsfinite remaining (not-yet-meshed)
        # edeges
        return done, params

    @process_text_tags(dim=2)
    def transfinite_surface_2D(self, done, params, dimtags, *args, **kwargs):
        arrangement = kwargs.get('arrangement', 'Left')
        cornerTags = kwargs.get('corner', [])

        for dim, tag in dimtags:
            if tag in done[2]:
                print("Surface " + str(tag) + " is already meshed (skipping)")
                continue
            gmsh.model.mesh.setTransfiniteSurface(tag,
                                                  arrangement=arrangement,
                                                  cornerTags=cornerTags)

        dimtags = [x for x in dimtags if not x[1] in done[2]]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(2)
        done[2].extend([x for dim, x in dimtags])
        return done, params

    @process_text_tags(dim=2)
    def transfinite_surface_3D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    transfinite_face_0D = transfinite_surface_0D
    transfinite_face_1D = transfinite_surface_1D
    transfinite_face_2D = transfinite_surface_2D
    transfinite_face_3D = transfinite_surface_3D

    # transfinite_volume
    @process_text_tags(dim=2)
    def transfinite_volume_0D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    @process_text_tags(dim=2)
    def transfinite_volume_1D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    @process_text_tags(dim=2)
    def transfinite_volume_2D(self, done, params, dimtags, *args, **kwargs):
        return done, params

    @process_text_tags(dim=2)
    def transfinite_volume_3D(self, done, params, dimtags, *args, **kwargs):
        arrangement = kwargs.get('arrangement', 'Left')
        cornerTags = kwargs.get('cornerTags', [])

        for dim, tag in dimtags:
            if tag in done[2]:
                print("Surface " + str(tag) + " is already meshed (skipping)")
                continue
            gmsh.model.mesh.setTransfiniteVolume(tag,
                                                 arrangement=arrangement,
                                                 cornerTags=cornerTags)

        dimtags = [x for x in dimtags if not x[1] in done[2]]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(3)
        done[3].extend([x for dim, x in dimtags])
        return done, params

    # merge text
    def _merge_xdim(self, mydim,  *args, **kwargs):
        text = kwargs.pop("text")
        dim = kwargs.pop("dim")
        if not dim[mydim]:
            return
        self.merge_text(text)

    def mergetxt_0D(self, done, params, *args, **kwargs):
        self._merge_xdim(0,  *args, **kwargs)
        return done, params

    def mergetxt_1D(self, done, params, *args, **kwargs):
        self._merge_xdim(1,  *args, **kwargs)
        return done, params

    def mergetxt_2D(self, done, params, *args, **kwargs):
        self._merge_xdim(2,  *args, **kwargs)
        return done, params

    def mergetxt_3D(self, done, params,  *args, **kwargs):
        self._merge_xdim(3,  *args, **kwargs)
        return done, params

    # copy edge
    @process_text_tags_sd(dim=2)
    def copy_edge_0D(self, vtag, tag1, tag2, nlayers):
        pass

    @process_text_tags_sd(dim=2)
    def copy_edge_1D(self, vtag, tag1, tag2, nlayers):
        pass

    @process_text_tags_sd(dim=2)
    def copy_edge_2D(self, vtag, tag1, tag2, nlayers):
        pass

    @process_text_tags_sd(dim=2)
    def copy_edge_3D(self, vtag, tag1, tag2, nlayers):
        pass

    # copy face
    @process_text_tags_sd(dim=2)
    def copyface_0D(self,  done, params, dimtags, dimtags2, *args, **kwargs):
        from petram.geom.geom_utils import find_translate_between_surface
        from petram.geom.geom_utils import find_rotation_between_surface
        from petram.geom.geom_utils import find_rotation_between_surface2

        ptx, p, l, s, v, mid_points = self.geom_info
        geom_size = np.sqrt(np.sum((np.max(ptx[:,0], 0) -  np.min(ptx[:,0], 0))**2))
        
        axan = kwargs.pop('axan', None)
        revolve = kwargs.pop('revolve', False)
        volume_hint = kwargs.pop('volume_hint', None)
        copy_cl = kwargs.pop('copy_cl', True)
        #geom_data = (ptx, l, s, v, mid_points)
        tag1 = [x for dim, x in dimtags]
        tag2 = [x for dim, x in dimtags2]

        if revolve:
            if volume_hint is None:
                # find volume hint
                # if there are volumes connecting all src and dst. we use thse volumes
                # to constuct point mapping
                volume_hint = []
                for v1 in v:
                    if (len(set(v[v1]).intersection(tag1)) != 0 and
                            len(set(v[v1]).intersection(tag2)) != 0):
                        volume_hint.append(str(v1))
                if len(volume_hint) == len(tag1):
                    volume_hint = ','.join(volume_hint)

            if volume_hint is None:
                ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = find_rotation_between_surface(
                    tag1, tag2, self.edge_tss,
                    geom_data=self.geom_info,
                    axan=axan, mind_eps=geom_size*self.mapper_tol)
            else:
                # copy(revolve) face is perfomed in the preparation of revolve mesh
                # in this case we use the volume being meshed as a hint
                #print("using volume hint", volume_hint)
                vtags = [int(x) for x in volume_hint.split(',')]
                #print("mind_eps", geom_size*self.mapper_tol)
                ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = find_rotation_between_surface2(
                    tag1, tag2, vtags, self.edge_tss,
                    geom_data=self.geom_info,
                    axan=axan, mind_eps=geom_size*self.mapper_tol)
        else:
            ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = find_translate_between_surface(
                tag1, tag2, self.edge_tss,
                geom_data=self.geom_info,
                axan=axan, mind_eps=geom_size*self.mapper_tol)

        params = (ax, an, px, d, affine, p_pairs, l_pairs, s_pairs)

        #print("p_pairs", p_pairs)
        for p0 in p_pairs:
            p1 = p_pairs[p0]
            if (0, p0) in self.cl_data:
                self.setCL(((0, p1),), self.cl_data[(0, p0)])
                #print("copying CL from ", p0, "to ", p1, self.cl_data[(0, p0)])
        #print("transformation param", params)
        return done, params

    @process_text_tags_sd(dim=2)
    def copyface_1D(self,  done, params, dimtags, dimtags2, *args, **kwargs):

        # don't do this otherwise previously meshed surfaces would be lost
        # gmsh.model.occ.synchronize()
        gmsh.model.mesh.rebuildNodeCache()
        ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = params

        #print("Entering CopyFace1D", l_pairs)

        ents = list(set(gmsh.model.getBoundary(
            dimtags, combined=False, oriented=False)))
        mdata = get_nodes_elements(ents, normalize=True)

        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)

        R = affine[:3, :3]
        D = affine[:3, -1]

        info1 = self.geom_info

        node_map2 = {}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info1[1][p_pairs[p]]+1

        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 1:
                continue

            tag = abs(l_pairs[tag])

            if tag in done[1]:
                print("Line is already meshed (CopyFace1D skip edge): " +
                      str(tag) + "... continuing")
                continue

            ntag, pos, ppos = ndata
            #print("parametric copyied", ppos)
            etypes, etags, nodes = edata
            ntag2 = range(noffset, noffset+len(ntag)-2)
            noffset = noffset+len(ntag)-2
            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])

            pos = np.array(pos).reshape(-1, 3).transpose()
            pos = (np.dot(R, pos).transpose() + D)
            pos = pos[:-2, :].flatten()

            # if node_map2[nodes[0][0]] == ptags2[0]:
            #    pass
            if len(nodes) > 1:
                assert False, "Only support linear geometry"

            # map start and end and check its parametricCoords
            nodes2 = [node_map2[nodes[0][0]], node_map2[nodes[0][-1]]]

            tmp = gmsh.model.mesh.getNodes(
                1, tag, includeBoundary=True)[2][-2:]
            p1_0 = gmsh.model.getValue(1, tag, [tmp[0]])
            p1_1 = gmsh.model.getValue(1, tag, [tmp[1]])

            p2_0 = info1[0][info1[1][nodes2[0]]]
            p2_1 = info1[0][info1[1][nodes2[1]]]

            def get_dist(p1, p2):
                d = np.array(p1) - np.array(p2)
                return np.sum(d**2)
            '''
            if get_dist(p1_0, p2_0) > get_dist(p1_0, p2_1):
                print("reversing Coords for tag :",
                      tag, p1_0, p1_1, p2_0, p2_1)
                pos = np.array(pos).reshape(-1, 3)
                pos = np.flip(pos, 0)
                pos = pos.flatten()
                print(pos, ppos, tmp)

                #print("parametric fixed", pos, ppos)
                #ntag2 = list(reversed(ntag2))
                #ppos = np.array([abs(1-x) for x in ppos])                
                #ppos = (tmp[0]-tmp[1])*ppos + tmp[1]                
            #else:
            '''
            ppos = (tmp[1]-tmp[0])*ppos + tmp[0]

            for i, j in zip(ntag, ntag2):
                node_map2[i] = j
            nodes2 = [[node_map2[x] for x in item] for item in nodes]

            gmsh.model.mesh.addNodes(dim, tag, ntag2, pos, ppos)
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)

            done[1].append(tag)

        return done, params

    @process_text_tags_sd(dim=2)
    def copyface_2D(self,  done, params, dimtags, dimtags2, *args, **kwargs):
        ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = params
 
        #print("Entering CopyFace2D", s_pairs)
        
        mdata = []
        for dim, tag in dimtags:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))

        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)

        R = affine[:3, :3]
        D = affine[:3, -1]

        info1 = self.geom_info
        node_map2 = {}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info1[1][p_pairs[p]]+1

        # add edge mapping
        ents_1D = self.expand_dimtags(dimtags, return_dim=1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([list(x[0]) for x in tmp], [])
        pos_s = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)
        pos_s = (np.dot(R, pos_s.transpose()).transpose() + D)

        ents_1D = self.expand_dimtags(dimtags2, return_dim=1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([list(x[0]) for x in tmp], [])
        pos_d = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        tree = cKDTree(pos_d)
        void, idx = tree.query(pos_s)
        #idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        for i, nt in zip(idx, ntags_s):
            node_map2[nt] = ntags_d[i]

        for kk, d in enumerate(mdata):
            dim, tag, ndata, edata = d

            if dim != 2:
                continue
            tag = s_pairs[tag]

            if tag in done[2]:
                print("Face is already meshed (CopyFace): " +
                      str(tag) + "... continuing")
                continue

            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2):
                node_map2[i] = j

            pos = np.array(pos).reshape(-1, 3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()

            #gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, [])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]

            gmsh.model.mesh.addNodes(dim, tag, ntag2, pos, [])
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)
            done[2].append(tag)
        # gmsh.model.mesh.reclassifyNodes()
        return done, params

    @process_text_tags_sd(dim=2)
    def copyface_3D(self,  done, params, dimtags, dimgtag2, *args, **kwargs):
        return done, params

    def _extract_info_for_volume(self, vtags):
        ptx, p, l, s, v, m = self.geom_info
        vv = {k: v[k] for k in vtags}
        ss = {k: s[k] for k in set(sum(vv.values(), []))}
        ll = {k: l[k] for k in set(sum(ss.values(), []))}
        pp = {k: p[k] for k in set(sum(ll.values(), []))}
        mm = {k: m[k] for k in ll}
        return ptx, pp, ll, ss, vv, mm

    # extrudeface
    @process_text_tags_vsd(dim=3)
    def extrude_face_0D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        from petram.geom.geom_utils import find_translate_between_surface
        from petram.geom.geom_utils import find_rotation_between_surface2
        from petram.geom.geom_utils import map_points_in_geom_info
        from petram.geom.geom_utils import map_lines_in_geom_info
        from petram.geom.geom_utils import map_surfaces_in_geom_info
        from petram.geom.geom_utils import map_volumes_in_geom_info

        revolve = kwargs.pop('revolve', False)
        nlayers = kwargs.get('nlayers', 5)
        use_recombine = kwargs.get('use_recombine', False)
        axan = kwargs.pop('axan', None)

        ptx, p, l, s, v, mid_points = self.geom_info
        geom_size = np.sqrt(np.sum((np.max(ptx[:,0], 0) -  np.min(ptx[:,0], 0))**2))

        #geom_data = (ptx, l, s, v, mid_points)
        tag1 = [x for dim, x in dimtags]
        tag2 = [x for dim, x in dimtags2]
        vtags = [x for dim, x in vdimtags]

        info1 = self._extract_info_for_volume(vtags)

        if revolve:
            ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = find_rotation_between_surface2(
                tag1, tag2, vtags, self.edge_tss,
                geom_data=self.geom_info,
                axan=axan, mind_eps=geom_size*self.mapper_tol)
        else:
            ax, an, px, d, affine, p_pairs, l_pairs, s_pairs = find_translate_between_surface(
                tag1, tag2, self.edge_tss,
                geom_data=self.geom_info,
                axan=axan, mind_eps=geom_size*self.mapper_tol)
        ws = self.prep_workspace()

        # ??? Apparently I need to delete a destinaion volume to make space for
        #     extrusion. Doesnt need to delete everything ???
        # self.delete_all_except_face(tag1)
        gmsh.model.occ.remove(vdimtags, recursive=False)

        if (not revolve and an == 0):
            ret = gmsh.model.occ.extrude(dimtags, d[0], d[1], d[2],
                                         numElements=[nlayers],
                                         recombine=use_recombine)
        elif (revolve and an != 0):
            ret = gmsh.model.occ.revolve(dimtags, px[0], px[1], px[2],
                                         ax[0], ax[1], ax[2], an,
                                         numElements=[nlayers],
                                         recombine=use_recombine)
        else:
            print(revolve, an)
            assert False, "extrude/revolve mesh error. Inconsistent imput"

        gmsh.model.occ.synchronize()

        # split dimtags to src, dst, lateral
        idx3 = np.where(np.array([x[0] for x in ret]) == 3)[0]
        vol = [ret[x] for x in idx3]
        dst = [ret[x-1] for x in idx3]
        laterals = [x for x in ret if not x in dst and x[0] == 2]
        ex_info = ([(2, x) for x in tag1], dst, laterals,  vol)

        # for debug the intermediate geometry
        #gmsh.write('tmp_'+ws +  '.brep')

        #info1 = self.geom_info
        info2 = self.read_geom_info(ret)

        pmap, pmap_r = map_points_in_geom_info(info1, info2,
                                               th=geom_size*self.mapper_tol)
        lmap, lmap_r = map_lines_in_geom_info(info1, info2, pmap_r,
                                              th=geom_size*self.mapper_tol)
        smap, smap_r = map_surfaces_in_geom_info(info1, info2, lmap_r)
        vmap, vmap_r = map_volumes_in_geom_info(info1, info2, smap_r)

        maps = (pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r)

        self.show_all()
        gmsh.model.mesh.generate(0)

        params = (ax, an, px, d, affine, p_pairs, l_pairs,
                  maps, ws, info1, info2, ex_info)
        self.switch_model('main1')

        return done, params

    @process_text_tags_vsd(dim=3)
    def extrude_face_1D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh lateral edges. copy destination edges
        # and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        ents1 = list(set(gmsh.model.getBoundary(
            dimtags, combined=False, oriented=False)))
        ents2 = list(set(gmsh.model.getBoundary(
            dimtags2, combined=False, oriented=False)))

        all_edges = list(lmap)
        meshed_edges = [k for k in all_edges if k in done[1]]
        tobe_meshed_edges = [k for k in all_edges if not k in done[1]]

        mdata = []
        for tag in meshed_edges:
            ndata = gmsh.model.mesh.getNodes(1, tag)
            edata = gmsh.model.mesh.getElements(1, tag)
            mdata.append((1, tag, ndata, edata))

        self.switch_model(ws)

        lateral_edge_digtags = [(1, lmap[key]) for key in tobe_meshed_edges]
        self.show_only(lateral_edge_digtags)
        
        gmsh.model.mesh.generate(1)
        #gmsh.write("debug_1d.msh")
        
        node_map1 = {info1[1][k]+1: info2[1][pmap[k]]+1 for k in pmap}
        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)

        # copy 1D elements on the source and dest surface
        copied_line = []
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim == 2:
                break
            tag = abs(lmap[tag])
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2):
                node_map1[i] = j

            vtags = [x for xx, x in gmsh.model.getBoundary(((1, tag),))]

            etypes, etags, nodes = edata

            etags2 = [list(range(eoffset, eoffset+len(etags[0]))),]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]

            # do I need to pay attetion the direction on parametricCoords here???
            gmsh.model.mesh.addNodes(dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)
            copied_line.append(tag)

        # gather data to transfer
        ents = [(1, x) for x in list(lmap_r)]
        ents = [x for x in ents if not x[1] in copied_line]
        mdata = get_nodes_elements(ents, normalize=True)
        '''
        mdata = []
        for dim, tag in ents:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))
        '''

        # send back 1D data
        self.switch_model('main1')
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)
        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        node_map3 = {info1[1][k]+1: info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]: x for x in node_map3}

        for d in mdata:
            dim, tag0, ndata, edata = d
            if dim != 1:
                continue

            tag = abs(lmap_r[tag0])
            #print("processing :", dim, tag0, " --> ", tag)
            if tag in done[1]:
                continue

            ntag, pos, ppos = ndata
            etypes, etags, nodes = edata

            ntag2 = range(noffset, noffset+len(ntag)-2)
            noffset = noffset+len(ntag)-2
            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])

            # omit edges
            pos = np.array(pos).reshape(-1, 3)
            pos = pos[:-2, :].flatten()

            nodes2 = [node_map3[nodes[0][0]], node_map3[nodes[0][-1]]]

            tmp = gmsh.model.mesh.getNodes(
                1, tag, includeBoundary=True)[2][-2:]
            p1_0 = gmsh.model.getValue(1, tag, [tmp[0]])
            p1_1 = gmsh.model.getValue(1, tag, [tmp[1]])
            p2_0 = info1[0][info1[1][nodes2[0]]]
            p2_1 = info1[0][info1[1][nodes2[1]]]

            def get_dist(p1, p2):
                d = np.array(p1) - np.array(p2)
                return np.sum(d**2)
            if get_dist(p1_0, p2_0) > get_dist(p1_0, p2_1):
                print("fixing parametricCoords for tag :",
                      tag, p1_0, p1_1, p2_0, p2_1)
                ppos = np.array([abs(1-x) for x in ppos])
                #ntag2 = list(reversed(ntag2))

            ppos = (tmp[1]-tmp[0])*ppos + tmp[0]

            for i, j in zip(ntag, ntag2):
                node_map3[i] = j
            nodes2 = [[node_map3[x] for x in item] for item in nodes]
            gmsh.model.mesh.addNodes(dim, tag, ntag2, pos, ppos)
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)

            done[1].append(tag)

        return done, params

    @process_text_tags_vsd(dim=3)
    def extrude_face_2D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh lateral face. copy destination face
        # and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        all_faces = list(smap)
        meshed_faces = [k for k in all_faces if k in done[2]]
        tobe_meshed_faces = [k for k in all_faces if not k in done[2]]

        # gather source/dest info
        mdata = []
        for tag in meshed_faces:
            ndata = gmsh.model.mesh.getNodes(2, tag)
            edata = gmsh.model.mesh.getElements(2, tag)
            mdata.append((2, tag, ndata, edata))

        # add edge mapping from main1
        ents_1D = self.expand_dimtags(vdimtags, return_dim=1)
        #ents_1D = self.expand_dimtags(list(dimtags)+list(dimtags2), return_dim = 1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([list(x[0]) for x in tmp], [])
        pos_s = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        self.switch_model(ws)
        src_dst = ex_info[0] + ex_info[1]
        src_dst_side = ex_info[0] + ex_info[1]+ex_info[2]
        # self.show_all()
        # self.hide(src_dst)
        self.show_only([(2, smap[key]) for key in tobe_meshed_faces])
        # this meshes un-meshed  sides ....

        gmsh.model.mesh.generate(2)
        
        ents_1D = gmsh.model.getEntities(1)
        #ents_1D = list(set(gmsh.model.getBoundary(src_dst, combined=False, oriented=False)))
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([list(x[0]) for x in tmp], [])
        pos_d = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        node_map1 = {info1[1][k]+1: info2[1][pmap[k]]+1 for k in pmap}

        #idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        tree = cKDTree(pos_d)
        void, idx = tree.query(pos_s)

        for i, nt in zip(idx, ntags_s):
            node_map1[nt] = ntags_d[i]
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)
        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # copy 2D elements on the source/destination surface
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 2:
                break
            tag = smap[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2):
                node_map1[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]

            gmsh.model.mesh.addNodes(dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)

        # self.show_all()
        # gmsh.write("debug_2d.msh")

        gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)
        self.show_only([(3, x) for x in list(vmap_r)])
        gmsh.model.mesh.generate(3)
        # gmsh.write("debug_3d.msh")

        # collect 3D mesh data
        tmp = [(3, x) for x in list(vmap_r)]
        mdata3D = []
        for dim, tag in tmp:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata3D.append((dim, tag, ndata, edata))

        # copy back lateral mesh on surfaces
        tmp = [(2, x) for x in list(smap_r)]
        tmp = [x for x in tmp if not x in src_dst]
        #print("gather info from ", tmp)
        mdata = []
        for dim, tag in tmp:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))

        ents_1D = gmsh.model.getEntities(1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([list(x[0]) for x in tmp], [])
        pos_s = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        # gmsh.write("debug.msh")
        self.switch_model('main1')

        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)

        ents_1D = self.expand_dimtags(vdimtags, return_dim=1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([list(x[0]) for x in tmp], [])
        pos_d = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        node_map3 = {info1[1][k]+1: info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]: x for x in node_map3}

        #idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        tree = cKDTree(pos_d)
        void, idx = tree.query(pos_s)

        for i, nt in zip(idx, ntags_s):
            node_map3[nt] = ntags_d[i]

        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 2:
                break
            tag = smap_r[tag]
            if tag in done[2]:
                continue

            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2):
                node_map3[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map3[x] for x in item] for item in nodes]

            gmsh.model.mesh.addNodes(
                dim, tag, ntag2, ndata[1], [])  # ndata[2])
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)
            done[2].append(tag)

        params = ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info, mdata3D
        return done, params

    @process_text_tags_vsd(dim=3)
    def extrude_face_3D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh volume and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info, mdata = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        # return done, params

        self.switch_model(ws)

        ents_1D = list(gmsh.model.getEntities(1)) + \
            list(gmsh.model.getEntities(2))
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([list(x[0]) for x in tmp], [])
        pos_s = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        self.switch_model('main1')
        #noffset = max(gmsh.model.mesh.getNodes()[0])+1
        #eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        noffset = int(max(gmsh.model.mesh.getNodes()[0])+1)
        eoffset = int(max(np.hstack(gmsh.model.mesh.getElements()[1]))+1)

        ents_1D = self.expand_dimtags(
            vdimtags, return_dim=1)+self.expand_dimtags(vdimtags, return_dim=2)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([list(x[0]) for x in tmp], [])
        pos_d = np.array(sum([list(x[1]) for x in tmp], [])).reshape(-1, 3)

        node_map3 = {info1[1][k]+1: info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]: x for x in node_map3}

        #idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        tree = cKDTree(pos_d)
        void, idx = tree.query(pos_s)

        for i, nt in zip(idx, ntags_s):
            node_map3[nt] = ntags_d[i]

        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 3:
                break
            tag = vmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2):
                node_map3[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map3[x] for x in item] for item in nodes]

            gmsh.model.mesh.addNodes(dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.addElements(dim, tag, etypes, etags2, nodes2)
            done[3].append(tag)
        return done, params

    # revolve
    # internally it is the same as extrude face
    def revolve_face_0D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve'] = True
        return self.extrude_face_0D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)

    def revolve_face_1D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve'] = True
        return self.extrude_face_1D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)

    def revolve_face_2D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve'] = True
        return self.extrude_face_2D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)

    def revolve_face_3D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve'] = True
        return self.extrude_face_3D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)

class GMSHMeshGeneratorBase():
    def __init__(self, q, task_q):
        self.q = q
        self.task_q = task_q

    def run(self):
        while True:
            time.sleep(0.1)
            try:
                task = self.task_q.get(True)
                self.ready_for_next_task()
            except EOFError:
                self.result_queue.put((-1, None))
                # self.task_queue.task_done()
                continue

            if task[0] == -1:
                # self.task_queue.task_done()
                break
            if task[0] == 1:
                try:
                    self.generate_mesh(*task[1])
                except BaseException:
                    txt = traceback.format_exc()
                    traceback.print_exc()
                    self.q.put((True, ('fail', txt)))
                    # self.task_queue.task_done()
                    break
        print("exiting prcesss")

    def generate_mesh(self, brep_input, msh_file, sequence,
                      dim, finalize, kwargs):

        kwargs['queue'] = self.q
        kwargs['mesh_sequence'] = sequence
        mw = GMSHMeshWrapper(**kwargs)

        max_dim, done, msh_output = mw.generate(brep_input, msh_file,
                                                dim=dim, finalize=finalize)

        self.q.put((True, (max_dim, done, msh_output)))


class GMSHMeshGenerator(GMSHMeshGeneratorBase, mp.Process):
    def __init__(self):

        # data to child
        task_q = mp.Queue()

        # data from child
        q = mp.Queue()

        GMSHMeshGeneratorBase.__init__(self, q, task_q)
        mp.Process.__init__(self)
        dprint1("starting a process for meshing")
        
    def ready_for_next_task(self):
        pass

from threading import Thread
from queue import Queue

class GMSHMeshGeneratorTH(GMSHMeshGeneratorBase, Thread):
    def __init__(self):

        # data to child
        task_q = Queue()

        # data from child
        q = Queue()

        GMSHMeshGeneratorBase.__init__(self, q, task_q)
        Thread.__init__(self)
        dprint1("starting a thread for mesh")        

    def ready_for_next_task(self):
        self.task_q.task_done()


