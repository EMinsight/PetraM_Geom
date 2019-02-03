from __future__ import print_function
import numpy as np
import gmsh
import time
import tempfile
import multiprocessing as mp
from Queue import Empty as QueueEmpty

from collections import OrderedDict
Algorithm2D= OrderedDict((("MeshAdap", 1), ("Automatic", 2), ("Delaunay", 3),
                ("Frontal", 6), ("BAMG", 7), ("DelQuad", 8),
                ("default", 2)))
Algorithm3D= OrderedDict((("Delaunay",1), ("New Delaunay",2),
                              ("Frontal", 4), 
                              ("Frontal Hex", 6), ("MMG3D", 7),
                              ("R-tree", 9), ("default", 1)))

import petram.geom.gmsh_config as gmsh_config

debug  = True
debug2 = False

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

def process_text_tags(dim=1, check = True):
    '''
    convert text tags input to dimtags
    tags = '1,2,3', 'all', or 'remaining'
    '''
    def func2(method):
        def method2(self, done, params, tags, *args, **kwargs):
            if tags == 'remaining':
                dimtags = gmsh.model.getEntities(dim)
                if check:
                    dimtags = [(dim, x) for xx, x in dimtags if not x in done[dim]]
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
    convert two text tags input to dimtags
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
        
class GMSHMeshWrapper(object):
    workspace_base = 1
    gmsh_init = False    
    def __init__(self, format=2.2,
                       CharacteristicLengthMax = 1e20,
                       CharacteristicLengthMin = 0,
                       EdgeResolution = 3, 
                       MeshAlgorithm = "Automatic",
                       MeshAlgorithm3D = "Delaunay",
                       MaxThreads = [1,1,1,1],
                       **kwargs):
        
        self.queue = kwargs.pop("queue", None)
        
        gmsh.clear()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", format)
        gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)

        gmsh_init = True
        self.add_model('main1')
        
        self.mesh_sequence = []
        
        # default clmax
        self.clmax = CharacteristicLengthMax
        self.clmin = CharacteristicLengthMin
        self.res = EdgeResolution 
        self.algorithm = MeshAlgorithm
        self.algorithm3d = MeshAlgorithm3D
        self.maxthreads= MaxThreads    # general, 1D, 2D, 3D: defualt = 1,1,1,1
        

        self._new_brep = True        

        
    def add(self, name, *gids, **kwargs):
        '''
        add mesh command
        '''
        if name == 'extrude_face':
            self.mesh_sequence.append(['copyface', (gids[1], gids[2]), kwargs])
        elif name == 'revolve_face':
            self.mesh_sequence.append(['copyface', (gids[1], gids[2]), kwargs])
        else:
            pass
        self.mesh_sequence.append([name, gids, kwargs])
        
    def count_sequence(self):
        return len(self.mesh_sequence)
    
    def clear(self):
        ''''
        clear mesh sequence
        '''
        self.mesh_sequence = []
        
    def generate(self, dim=3, brep_input = '',
                 finalize=False, msh_file=''):
        '''
        generate mesh based on  meshing job sequence.
        brep must be loaed 
        '''
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
        gmsh.option.setNumber("Mesh.Algorithm3D",
                              Algorithm3D[self.algorithm3d])
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",
                              self.clmax)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",
                              self.clmin)
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        gmsh.option.setNumber("General.NumThreads",   self.maxthreads[0])        
        gmsh.option.setNumber("Mesh.MaxNumThreads1D", self.maxthreads[1])
        gmsh.option.setNumber("Mesh.MaxNumThreads2D", self.maxthreads[2])
        gmsh.option.setNumber("Mesh.MaxNumThreads3D", self.maxthreads[3])

        # 
        self.vertex_geom_size = get_vertex_geom_zie()
        # set default vertex mesh size
        for tag in self.vertex_geom_size.keys():
            size = self.vertex_geom_size[tag]/self.res
            if size > self.clmax: size = self.clmax
            if size <= self.clmin: size = self.clmin
            gmsh.model.mesh.setSize(((0, tag),), size)
        
        done = [[], [], [], []]
        params = [None]*len(self.mesh_sequence)

        maxdim = max([x for x, tag in gmsh.model.getEntities()])
        maxdim = min([maxdim, dim])
        for mdim in range(maxdim+1):
            for idx, sq in enumerate(self.mesh_sequence):
                proc, args, kwargs = sq
                f = getattr(self, proc+"_"+str(mdim)+"D")
                
                if self.queue is not None:
                    self.queue.put((False,
                                    "Processing " + proc+"_"+str(mdim)+"D"))
                
                done, params[idx] = f(done, params[idx],
                                      *args, **kwargs)
                
            for i in range(mdim+1, 4): done[i] = []
            
        if finalize:
            self.add_sequential_physicals()
            if msh_file != '': gmsh.write(msh_file)

        return maxdim, done
    
    def load_brep(self, filename):
        self.switch_model('main1')
        self.target = filename
        
        gmsh.model.occ.importShapes(filename, highestDimOnly=False)
        gmsh.model.occ.synchronize()

        self.geom_info = self.read_geom_info()
        self._new_brep = True
        return self.current
    
    def read_geom_info(self):
        from petram.geom.read_gmsh import read_loops2
        self.hide_all()        
        gmsh.model.mesh.generate(1)

        #ptx, p, l, s, v = read_loops2(gmsh)
        return read_loops2(gmsh)
        
    def add_model(self, name):
        gmsh.model.add(name)
        name = gmsh.model.list()[-1]
        self.current = name
        return name
        
    def switch_model(self, name = 'main1'):
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
        gmsh.model.setVisibility(dimtags, True, recursive = recursive)
        
    def hide(self, dimtags, recursive=False):
        gmsh.model.setVisibility(dimtags, False, recursive = recursive)
        
    def hide_all(self):
        ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, False)
        
    def show_all(self):
        ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, True)

    def delete_all_except(self, dim, tags):
        for d in [3, 2, 1]:
            ent = gmsh.model.getEntities(d)            
            if dim == d:
                ent = [(d, t) for d, t in ent if not t in tags]
                gmsh.model.occ.remove(ent, recursive=True)
                break
            else:
                gmsh.model.occ.remove(ent)
        gmsh.model.occ.synchronize()        
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

    def expand_dimtags(self, dimtags, return_dim = None):
        '''
        expand all subdims
        '''
        volumes = []
        surfaces = []
        lines = []
        points = []
        
        for dim, tag in dimtags:
            if dim == 3: volumes.append(tag)
            if dim == 2: surfaces.append(tag)
            if dim == 1: lines.append(tag)
            if dim == 0: points.append(tag)
            
        for v in volumes:           
            surfaces.extend(self.geom_info[4][v])

        surfaces = list(set(surfaces))
        for s in surfaces:
            lines.extend(self.geom_info[3][s])
                        
        lines = list(set(lines))
        for l in lines:
            points.extend(self.geom_info[2][l])

        points = list(set(points))

        ans = ( [(0, p) for p in points] +
                [(1, l) for l in lines] +
                [(2, s) for s in surfaces] +
                [(3, v) for v in volumes] )
        
        if return_dim is not None:
            return [(return_dim, x[1]) for x in ans if x[0]==return_dim]
        
        return ans

    def list_entities(self):
        for x in gmsh.model.list():
            gmsh.model.setCurrent(x)
            print(x, gmsh.model.getEntities())
        gmsh.model.setCurrent(self.current)

    def add_sequential_physicals(self, verbose = True, include_lower_dims=False):
        '''
        add sequencial physical entity numbers
        '''
        if self.queue is not None:
            self.queue.put((False, "Adding Physicals..."))
            
        ent = gmsh.model.getEntities(dim=3)
        max_dim = 0
        
        if len(ent) > 0:
            max_dim = 3
            if verbose: print("Adding " + str(len(ent)) + " Volume(s)")
        for k, x in enumerate(ent):
            if len(gmsh.model.getPhysicalGroupsForEntity(3, x[1])) > 0:
                continue
            value = gmsh.model.addPhysicalGroup(3, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(3, value, 'volume'+str(value))
                  
        ent = gmsh.model.getEntities(dim=2)
        if len(ent) > 0:
            if max_dim == 0: max_dim = 2                              
            if verbose: print("Adding " + str(len(ent)) + " Surface(s)")
            
        for k, x in enumerate(ent):
            if len(gmsh.model.getPhysicalGroupsForEntity(2, x[1])) > 0:
                continue
            value = gmsh.model.addPhysicalGroup(2, [x[1]], tag=k+1)
            gmsh.model.setPhysicalName(2, value, 'surface'+str(value))

        if not include_lower_dims and max_dim == 3: return
        ent = gmsh.model.getEntities(dim=1)
        if len(ent) > 0:
            if max_dim == 0: max_dim = 1
            if verbose: print("Adding " + str(len(ent)) + " Line(s)")
            
        for k, x in enumerate(ent):
            if len(gmsh.model.getPhysicalGroupsForEntity(1, x[1])) > 0:
                continue                    
            value = gmsh.model.addPhysicalGroup(1, [x[1]], tag=k+1)                
            gmsh.model.setPhysicalName(1, value, 'line'+str(value))
            
        if not include_lower_dims and max_dim == 2: return                      
        ent = gmsh.model.getEntities(dim=0)
        if len(ent) > 0:
            if verbose: print("Adding " + str(len(ent)) + " Point(s)")
        for k, x in enumerate(ent):
            if len(gmsh.model.getPhysicalGroupsForEntity(0, x[1])) > 0:
                continue                                        
            value = gmsh.model.addPhysicalGroup(0, [x[1]], tag=k+1)                
            gmsh.model.setPhysicalName(0, value, 'point'+str(value))
        
    '''
    Low-level implementation at each mesh dim
    '''
    @process_text_tags(dim=0)            
    def cl_0D(self, done, params, dimtags, size, overwrite = True):
        '''
        Set Charcteristic length of points
        By default, it overwrite whatever it is set
        '''
        
        if not overwrite:
            dimtags = [(0, x) for x in tags]        
        else:
            dimtags = [x for x in dimtags if not x[0] in done[0]]
        if len(dimtags) == 0: return
        self.show_only(dimtags)
        gmsh.model.mesh.setSize(dimtags, size)
        for dim, tag in dimtags:
            if not tag in done[0]:  done[0].append(tag)
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
    def recombine_face_0D(self, done, params, dimtags):
        for dim, tag in dimtags:
            if dim != 2: continue
            gmsh.model.mesh.setRecombine(2, tag)
        return done, params
    
    @process_text_tags(dim=2, check=False)
    def recombine_face_1D(self, done, params, dimtags):
        return done, params
    
    @process_text_tags(dim=2, check=False)
    def recombine_face_2D(self, done, params, dimtags):
        return done, params
    
    @process_text_tags(dim=2, check=False)
    def recombine_face_3D(self, done, params, dimtags):
        return done, params        

    # freevolume
    @process_text_tags(dim=3)        
    def freevolume_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize=kwargs.pop("maxsize", 1e20)
        minsize=kwargs.pop("minsize", 0.0)
        res=kwargs.pop("resolution",  np.inf)

        embeds = [x for x in kwargs.pop("embeds",  '').split(',')]
        embeds = [x for x in embeds if len(x) > 0]                        
        if len(embeds) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embeds:
                gmsh.mesh.embed(2, x, dimtags[0][0], dimtags[0][1])
        embedl=[x for x in kwargs.pop("embedl",  '').split(',')]
        embedl = [x for x in embedl if len(x) > 0]                
        if len(embedl) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embedl:
                gmsh.mesh.embed(2, x, dimtags[0][0], dimtags[0][1])
        embedp = [x for x in kwargs.pop("embedp",  '').split(',')]
        embedp = [x for x in embedp if len(x) > 0]                                
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embedp:
                gmsh.mesh.embed(1, x, dimtags[0][0], dimtags[0][1])

        dimtags = self.expand_dimtags(dimtags, return_dim = 0)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize: size = maxsize
            if size < minsize: size = minsize
            gmsh.model.mesh.setSize(((0, tag),), size)
            done[0].append(tag)            
        gmsh.model.mesh.generate(0)
        return done, params
    
    @process_text_tags(dim=3)    
    def freevolume_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[3].extend([x for dim, x in dimtags])
        sdimtags = self.expand_dimtags(dimtags, return_dim = 2)
        done[2].extend([x for dim, x in sdimtags if not x in done[2]])        
        
        dimtags = self.expand_dimtags(dimtags, return_dim = 1)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]
        tags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]        
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in tags])        
        return done, params
    
    @process_text_tags(dim=3)    
    def freevolume_2D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[3].extend([x for dim, x in dimtags])
        
        dimtags = self.expand_dimtags(dimtags, return_dim = 2)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[2]]
        tags = [(dim, tag) for dim, tag in dimtags if not tag in done[2]]

        self.show_only(dimtags)
        gmsh.model.mesh.generate(2)
        done[2].extend([x for dim, x in tags])                
        return done, params
    
    @process_text_tags(dim=3)
    def freevolume_3D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(3)
        done[3].extend([x for dim, x in tags])                        
        return done, params
    
    # freeface
    @process_text_tags(dim=2)
    def freeface_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize=kwargs.pop("maxsize", 1e20)
        minsize=kwargs.pop("minsize", 0.0)
        res=kwargs.pop("resolution", np.inf)

        embedl=[x for x in kwargs.pop("embedl",  '').split(',')]
        embedl = [x for x in embedl if len(x) > 0]        
        if len(embedl) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embedl:
                gmsh.mesh.embed(2, x, dimtags[0][0], dimtags[0][1])
        embedp=[x for x in kwargs.pop("embedp",  '').split(',')]
        embedp = [x for x in embedp if len(x) > 0]        
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embedp:
                gmsh.mesh.embed(1, x, dimtags[0][0], dimtags[0][1])
                
        dimtags = self.expand_dimtags(dimtags, return_dim = 0)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize: size = maxsize
            if size < minsize: size = minsize
            gmsh.model.mesh.setSize(((0, tag),), size)
            done[0].append(tag)            
        gmsh.model.mesh.generate(0)
        return done, params
    
    @process_text_tags(dim=2)                
    def freeface_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)

        done[2].extend([x for dim, x in dimtags])
        
        dimtags = self.expand_dimtags(dimtags, return_dim = 1)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]
        tags = [(dim, tag) for dim, tag in dimtags if not tag in done[1]]        
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in tags])                                
        return done, params
    
    @process_text_tags(dim=2)                
    def freeface_2D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(2)
        done[2].extend([x for dim, x in tags])         
        return done, params

    @process_text_tags(dim=2)    
    def freeface_3D(self, done, params, tags, *args, **kwargs):
        return done, params                

    # freeedge
    @process_text_tags(dim=1)        
    def freeedge_0D(self, done, params, dimtags, *args, **kwargs):
        maxsize=kwargs.pop("maxsize", 1e20)
        minsize=kwargs.pop("minsize", 0.0)
        res=kwargs.pop("resolution", np.inf)

        embedp=[x for x in kwargs.pop("embedp",  '').split(',')]
        embedp = [x for x in embedp if len(x) > 0]
        if len(embedp) > 0:
            if len(dimtags) > 1:
                assert False, "Embed works only when there is one target"
            for x in embedp:
                gmsh.mesh.embed(1, x, dimtags[0][0], dimtags[0][1])
                
        done[1].extend([x for dim, x in dimtags])
        
        dimtags = self.expand_dimtags(dimtags, return_dim = 0)
        dimtags = [(dim, tag) for dim, tag in dimtags if not tag in done[0]]
        self.show_only(dimtags)
        for dim, tag in dimtags:
            size = self.vertex_geom_size[tag]/res
            if size > maxsize: size = maxsize
            if size <= minsize: size = minsize
            gmsh.model.mesh.setSize(((0, tag),), size)
            done[0].append(tag)            
        gmsh.model.mesh.generate(0)
        return done, params
    
    @process_text_tags(dim=1)                
    def freeedge_1D(self, done, params, dimtags, *args, **kwargs):
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        tags = [(dim, tag) for dim, tag in dimtags]
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in tags])                 
        return done, params

    @process_text_tags(dim=1)    
    def freeedge_2D(self, done, params, tags, *args, **kwargs):
        return done, params        

    @process_text_tags(dim=1)    
    def freeedge_3D(self,  done, params, tags, *args, **kwargs):
        return done, params                

    # transfinite_edge
    @process_text_tags(dim=1)        
    def transfinite_edge_0D(self, done, params, dimtags, *args, **kwargs):
        #meher.add('transfinite_line', gid, nseg=nseg, progression = p,  bump = b)
        nseg = kwargs.get('nseg', 100)-1
        
        meshType = 'Progression'
        coef = kwargs.get('progression', 1.0)
        
        bump = kwargs.get('bump', 1)
        if bump  != 1:
            meshType = 'Bump'
            coef = bump

        for dim, tag in dimtags:
            if tag in done[1]: continue
            gmsh.model.mesh.setTransfiniteCurve(tag, nseg, 
                                                meshType = meshType,
                                                coef = coef)
            done[1].append(tag)
        return done, params

        '''
        # A test to understand transfinite meshing...
        if bump != 1:
            if nseg % 2 == 1:            
                tmp = np.logspace(0, np.log10(bump), (nseg+1)/2)
                dd =  np.hstack([np.flip(tmp[1:], 0), tmp])
            else:
                tmp = np.logspace(0, np.log10(bump), nseg/2)
                dd =  np.hstack([np.flip(tmp, 0), tmp])
        elif progression != 1:
             dd = abs(progression)**np.arange(nseg)
             if progression < 0:
                 dd = np.flip(dd, 0)
        else:
            dd = np.array([1.]*nseg)

        dists = dd/np.sum(dd) # normalized goemetrical distance of nodes
        cdists = np.hstack([[0], np.cumsum(dists)])
        

        params = {}
        for dim, tag in dimtags:
            ntags, nodes, pcoords = gmsh.model.mesh.getNodes(1, tag, includeBoundary=True)
            p = list(np.linspace(pcoords[-2], pcoords[-1], np.max([nseg*3, 300])))            
            ptx = np.array(gmsh.model.getValue(1, tag, p)).reshape(-1, 3)
            dd = np.sqrt(np.sum((ptx[:-1]- ptx[1:])**2, 1))
            ddd = np.hstack([[0], np.cumsum(dd)])
            ddd = ddd/ddd[-1]

            pcoords = np.interp(cdists, ddd, p) #parametric Coords to realized the distancs
            nodepos = gmsh.model.getValue(1, tag, pcoords[1:-1])

            size = np.array(gmsh.model.getValue(1, tag, pcoords[:2])).reshape(-1, 3)
            size1 = np.sqrt(np.sum((size[1]- size[0])**2))
            size = np.array(gmsh.model.getValue(1, tag, pcoords[-2:])).reshape(-1, 3)
            size2 = np.sqrt(np.sum((size[1]- size[0])**2))
                
            vtags = [x for xx, x in gmsh.model.getBoundary(((dim, tag),))]
            flag = check_line_orientation(tag, vtags, pcoords[0])            
            if flag == 1:
                if not vtags[0] in done[0]:                
                    gmsh.model.mesh.setSize(((0, vtags[0]),), size1)
                    done[0].append(vtags[0])                                
                if not vtags[1] in done[0]:                                    
                    gmsh.model.mesh.setSize(((0, vtags[1]),), size2)
                    done[0].append(vtags[1])
                    
            elif flag == 2:
                if not vtags[0] in done[0]:
                    gmsh.model.mesh.setSize(((0, vtags[0]),), size2)
                    done[0].append(vtags[0])                                                    
                if not vtags[1] in done[0]:                    
                    gmsh.model.mesh.setSize(((0, vtags[1]),), size1)
                    done[0].append(vtags[1])                    
            else:
                print(gmsh.model.getValue(0, vtags[0], [0]), gmsh.model.getValue(0, vtags[1], [0]),
                      gmsh.model.getValue(1, tag, pcoords[:1]), gmsh.model.getValue(1, tag, pcoords[-1:]))
                assert False, "Something is wrong"
            
            params[tag] = (nodepos, pcoords)
        gmsh.model.mesh.generate(0)
        return done, params
        '''
    @process_text_tags(dim=1)            
    def transfinite_edge_1D(self, done, params, dimtags, *args, **kwargs):
        self.show_only(dimtags)
        gmsh.model.mesh.generate(1)
        done[1].extend([x for dim, x in dimtags])         
        return done, params
        
        '''
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        
        for dim, tag in dimtags:
            nodepos, pcoords = params[tag]
            
            ntags = range(noffset, noffset+len(pcoords)-2)
            noffset = noffset+len(ntags)
            
            etags = [range(eoffset, eoffset+len(ntags)+1)]
            eoffset = eoffset+len(etags[0])

            vtags = [x for xx, x in gmsh.model.getBoundary(((dim, tag),))]
            flag = check_line_orientation(tag, vtags, pcoords[0])
            
            info1 = self.geom_info
            tmp = np.vstack([ntags, ntags]).transpose().flatten()
            if flag == 1:
                nodes2 = np.hstack([[info1[1][vtags[0]]+1], tmp, [info1[1][vtags[1]]+1]])
            elif flag == 2:
                nodes2 = np.hstack([[info1[1][vtags[1]]+1], tmp, [info1[1][vtags[0]]+1]])                
            else:
                assert False, "Something is wrong"

            #print("setting", dim, tag, ntags, nodepos, pcoords)
            gmsh.model.mesh.setNodes(dim, tag, ntags, nodepos, pcoords[1:-1])
            #print("setting", dim, tag, [1], etags, [list(nodes2)])     
            gmsh.model.mesh.setElements(dim, tag, [1], etags, [list(nodes2)])        
            done[1].append(tag)
            
        return done, params                        
        '''
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
                print("Surface "+ str(tag) + " is already meshed (skipping)")
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
                print("Surface "+ str(tag) + " is already meshed (skipping)")
                continue
            gmsh.model.mesh.setTransfiniteVolume(tag,
                                                  arrangement=arrangement,
                                                  cornerTags=cornerTags)
            
        dimtags = [x for x in dimtags if not x[1] in done[2]]                    
        self.show_only(dimtags)
        gmsh.model.mesh.generate(3)
        done[3].extend([x for dim, x in dimtags])
        return done, params        
        
    # copy edge    
    def copy_edge_0D(self, vtag, tag1, tag2, nlayers):
        pass
    def copy_edge_1D(self, vtag, tag1, tag2, nlayers):
        pass
    def copy_edge_2D(self, vtag, tag1, tag2, nlayers):
        pass
    def copy_edge_3D(self, vtag, tag1, tag2, nlayers):
        pass

    # copy face
    @process_text_tags_sd(dim=2)        
    def copyface_0D(self,  done, params, dimtags, dimtags2, *args, **kwargs):
        from petram.geom.geom_utils import find_translate_between_surface
        ptx, p, l, s, v = self.geom_info
        axan = kwargs.pop('axan', None)        
        geom_data = (ptx, l, s, None)
        tag1 = [x for dim, x in dimtags]
        tag2 = [x for dim, x in dimtags2]        
        ax, an, px, d, affine, p_pairs, l_pairs = find_translate_between_surface(
                                                            tag1, tag2,
                                                            geom_data=geom_data,
                                                            axan = axan)
        
        params = (ax, an, px, d, affine, p_pairs, l_pairs)
        return done, params
    
    @process_text_tags_sd(dim=2)
    def copyface_1D(self,  done, params, dimtags, dimtags2, *args, **kwargs):
        
        ax, an, px, d, affine, p_pairs, l_pairs = params
        
        ents = list(set(gmsh.model.getBoundary(dimtags, combined=False, oriented=False)))
        mdata = get_nodes_elements(ents, normalize=True)
        
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        R = affine[:3,:3]
        D = affine[:3,-1]

        info1 = self.geom_info        
        node_map2 = {}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info1[1][p_pairs[p]]+1

        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 1: continue

            tag = l_pairs[tag]
            if tag in done[1]:
                print("Line is already meshed (CopyFace failed): "+str(tag)+ "... continuing")
                continue
                
            ntag, pos, ppos = ndata            
            etypes, etags, nodes = edata
            ntag2 = range(noffset, noffset+len(ntag)-2)
            noffset = noffset+len(ntag)-2
            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            
            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D)
            pos = pos[:-2,:].flatten()            

            #if node_map2[nodes[0][0]] == ptags2[0]:
            #    pass
            if len(nodes)>1:
                assert False, "Only support linear geometry"

            # map start and end and check its parametricCoords
            nodes2 = [node_map2[nodes[0][0]], node_map2[nodes[0][-1]]]

            tmp = gmsh.model.mesh.getNodes(1, tag, includeBoundary=True)[2][-2:]
            p1_0 = gmsh.model.getValue(1, tag, [tmp[0]])
            p1_1 = gmsh.model.getValue(1, tag, [tmp[1]])            
            p2_0 = info1[0][info1[1][nodes2[0]]]
            p2_1 = info1[0][info1[1][nodes2[1]]]            

            def get_dist(p1, p2):
                d = np.array(p1) - np.array(p2)
                return np.sum(d**2)
            if get_dist(p1_0, p2_0) > get_dist(p1_0, p2_1):
                print("fixing parametricCoords for tag :", tag, p1_0, p1_1, p2_0, p2_1)
                ppos = np.array([abs(1-x) for x in ppos])
                #ntag2 = list(reversed(ntag2))

            ppos = (tmp[1]-tmp[0])*ppos + tmp[0]
          
            for i, j in zip(ntag, ntag2): node_map2[i] = j
            nodes2 = [[node_map2[x] for x in item] for item in nodes]
            
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ppos)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)        
            done[1].append(tag)
            
        return done, params            

    @process_text_tags_sd(dim=2)    
    def copyface_2D(self,  done, params, dimtags, dimtags2, *args, **kwargs):
        ax, an, px, d, affine, p_pairs, l_pairs = params
        
        mdata = []
        for dim, tag in dimtags:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))
            
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        R = affine[:3,:3]
        D = affine[:3,-1]

        info1 = self.geom_info        
        node_map2 = {}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info1[1][p_pairs[p]]+1

        ## add edge mapping    
        ents_1D = self.expand_dimtags(dimtags, return_dim = 1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([x[0] for x in tmp], [])
        pos_s = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        pos_s = (np.dot(R, pos_s.transpose()).transpose() + D)
        
        ents_1D = self.expand_dimtags(dimtags2, return_dim = 1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([x[0] for x in tmp], [])
        pos_d = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        
        idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        for i, nt in zip(idx, ntags_s):
            node_map2[nt] = ntags_d[i]
            

        for kk, d in enumerate(mdata):
            dim, tag, ndata, edata = d
            if dim != 2: continue
            tag = dimtags2[kk][1]
            if tag in done[2]:
                print("Face is already meshed (CopyFace): "+str(tag) + "... continuing")
                continue                
                
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map2[i] = j

            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()

            #gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, [])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]

            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)        
            done[2].append(tag)
            
        return done, params            
        

    @process_text_tags_sd(dim=2)        
    def copyface_3D(self,  done, params, dimtags, dimgtag2, *args, **kwargs):
        return done, params                            
    
    # extrudeface
    @process_text_tags_vsd(dim=3)            
    def extrude_face_0D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        from petram.geom.geom_utils import find_translate_between_surface
        from petram.geom.geom_utils import map_points_in_geom_info
        from petram.geom.geom_utils import map_lines_in_geom_info
        from petram.geom.geom_utils import map_surfaces_in_geom_info
        from petram.geom.geom_utils import map_volumes_in_geom_info
        
        revolve = kwargs.pop('revolve', False)
        nlayers = kwargs.get('nlayers', 5)
        axan = kwargs.pop('axan', None)
        
        ptx, p, l, s, v = self.geom_info
        geom_data = (ptx, l, s, None)
        tag1 = [x for dim, x in dimtags]
        tag2 = [x for dim, x in dimtags2]        
        ax, an, px, d, affine, p_pairs, l_pairs = find_translate_between_surface(
                                                            tag1, tag2,
                                                            geom_data=geom_data,
                                                            axan = axan)
        
        ws = self.prep_workspace()

        self.delete_all_except(2, tag1)
        if (not revolve and an == 0):
            ret = gmsh.model.occ.extrude(dimtags, d[0], d[1], d[2],
                                         numElements=[nlayers])
        elif (revolve and an != 0):
            ret = gmsh.model.occ.revolve(dimtags, px[0], px[1], px[2],
                                         ax[0], ax[1], ax[2], an, 
                                         numElements=[nlayers])
        else:
            print(revolve, an)
            assert False, "extrude/revolve mesh error. Inconsistent imput"

        # split dimtags to src, dst, lateral
        idx3 = np.where(np.array([x[0] for x in ret]) == 3)[0]
        vol = [ret[x] for x in idx3]
        dst = [ret[x-1] for x in idx3]
        laterals = [x for x in ret if not x in dst and x[0] == 2] 
        ex_info = ([(2, x) for x in tag1], dst, laterals,  vol)
        gmsh.model.occ.synchronize()
        
        # for debug the intermediate geometry
        gmsh.write('tmp_'+ws +  '.brep')            

        info1 = self.geom_info
        info2 = self.read_geom_info()
        gmsh.model.getEntities

        pmap, pmap_r = map_points_in_geom_info(info1, info2)
        lmap, lmap_r = map_lines_in_geom_info(info1, info2, pmap_r)
        smap, smap_r = map_surfaces_in_geom_info(info1, info2, lmap_r)
        vmap, vmap_r = map_volumes_in_geom_info(info1, info2, smap_r)

        maps = (pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r)
        self.show_all()
        gmsh.model.mesh.generate(0)

        params = (ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info)
        self.switch_model('main1')
        
        return done, params
    
    @process_text_tags_vsd(dim=3)                    
    def extrude_face_1D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh lateral edges. copy destination edges
        # and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        ents1 = list(set(gmsh.model.getBoundary(dimtags, combined=False, oriented=False)))
        ents2 = list(set(gmsh.model.getBoundary(dimtags2, combined = False, oriented=False)))
        mdata = []
        for dim, tag in ents1+ents2:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))

        self.switch_model(ws)

        self.show_all()
        gmsh.model.mesh.generate(1)
        
        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}        
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # copy 1D elements on the source and dest surface
        copied_line = []
        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim == 2: break
            tag = abs(lmap[tag])
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            vtags = [x for xx, x in gmsh.model.getBoundary(((1, tag),))]
            
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]

            # do I need to pay attetion the direction on parametricCoords here???

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])            
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            copied_line.append(tag)

        # gather data to transfer
        ents = gmsh.model.getEntities(1)
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
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        node_map3 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]:x for x in node_map3}

        for d in mdata:
            dim, tag0, ndata, edata = d
            if dim != 1: continue

            tag = abs(lmap_r[tag0])
            #print("copy from ", dim, tag0, "to", tag)
            
            ntag, pos, ppos = ndata
            etypes, etags, nodes = edata
            
            ntag2 = range(noffset, noffset+len(ntag)-2)
            noffset = noffset+len(ntag)-2
            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])

            # omit edges
            pos = np.array(pos).reshape(-1,3)
            pos = pos[:-2,:].flatten()            

            nodes2 = [node_map3[nodes[0][0]], node_map3[nodes[0][-1]]]
            
            tmp = gmsh.model.mesh.getNodes(1, tag, includeBoundary=True)[2][-2:]
            p1_0 = gmsh.model.getValue(1, tag, [tmp[0]])
            p1_1 = gmsh.model.getValue(1, tag, [tmp[1]])            
            p2_0 = info1[0][info1[1][nodes2[0]]]
            p2_1 = info1[0][info1[1][nodes2[1]]]            

            def get_dist(p1, p2):
                d = np.array(p1) - np.array(p2)
                return np.sum(d**2)
            if get_dist(p1_0, p2_0) > get_dist(p1_0, p2_1):
                print("fixing parametricCoords for tag :", tag, p1_0, p1_1, p2_0, p2_1)
                ppos = np.array([abs(1-x) for x in ppos])
                #ntag2 = list(reversed(ntag2))
                
            ppos = (tmp[1]-tmp[0])*ppos + tmp[0]
            
            for i, j in zip(ntag, ntag2): node_map3[i] = j                        
            nodes2 = [[node_map3[x] for x in item] for item in nodes]
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ppos)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            
            done[1].append(tag)
            
        return done, params
    
    @process_text_tags_vsd(dim=3)
    def extrude_face_2D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh lateral face. copy destination face
        # and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        # gather source/dest info
        mdata = []
        for dim, tag in list(dimtags)+list(dimtags2):
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))

        ## add edge mapping from main1   
        ents_1D = self.expand_dimtags(list(dimtags)+list(dimtags2), return_dim = 1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([x[0] for x in tmp], [])
        pos_s = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        
        self.switch_model(ws)
        src_dst = ex_info[0]+ ex_info[1]
        src_dst_side = ex_info[0]+ ex_info[1]+ex_info[2]        
        self.show_all()
        self.hide(src_dst)
        # this mesh sides....
        gmsh.model.mesh.generate(2)
        
        ents_1D = list(set(gmsh.model.getBoundary(src_dst, combined=False, oriented=False)))
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([x[0] for x in tmp], [])
        pos_d = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}                
        idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        for i, nt in zip(idx, ntags_s):
            node_map1[nt] = ntags_d[i]
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # copy 2D elements on the source/destination surface
        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim != 2: break
            tag = smap[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])            
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.show_all()
        gmsh.model.mesh.generate(3)

        tmp = gmsh.model.getEntities(3)
        mdata3D = []
        for dim, tag in tmp:
            ndata = gmsh.model.mesh.getNodes(dim, tag) 
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata3D.append((dim, tag, ndata, edata))

        # copy back lateral mesh on surfaces
        tmp = gmsh.model.getEntities(2)
        tmp = [x for x in tmp if not x in src_dst]
        #print("gather info from ", tmp)
        mdata = []
        for dim, tag in tmp:
            ndata = gmsh.model.mesh.getNodes(dim, tag) 
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))
        
        ents_1D = gmsh.model.getEntities(1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([x[0] for x in tmp], [])
        pos_s = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        self.switch_model('main1')
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        ents_1D = self.expand_dimtags(vdimtags, return_dim = 1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([x[0] for x in tmp], [])
        pos_d = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        
        node_map3 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]:x for x in node_map3}
        idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        for i, nt in zip(idx, ntags_s):
            node_map3[nt] = ntags_d[i]

        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim != 2: break
            tag = smap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map3[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map3[x] for x in item] for item in nodes]

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], [])#ndata[2])            
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            done[2].append(tag)
            
        params = ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info, mdata3D            
        return done, params        

    @process_text_tags_vsd(dim=3)
    def extrude_face_3D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        # mesh volume and bring them back to main1
        ax, an, px, d, affine, p_pairs, l_pairs, maps, ws, info1, info2, ex_info, mdata = params
        pmap, pmap_r, lmap, lmap_r, smap, smap_r, vmap, vmap_r = maps

        #return done, params
    
        self.switch_model(ws)

        ents_1D = list(gmsh.model.getEntities(1))+list(gmsh.model.getEntities(2))
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_s = sum([x[0] for x in tmp], [])
        pos_s = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        self.switch_model('main1')
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        ents_1D = self.expand_dimtags(vdimtags, return_dim = 1)+self.expand_dimtags(vdimtags, return_dim = 2)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_d = sum([x[0] for x in tmp], [])
        pos_d = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        
        node_map3 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map3 = {node_map3[x]:x for x in node_map3}
        idx = [np.argmin(np.sum((pos_d-p)**2, 1)) for p in pos_s]
        for i, nt in zip(idx, ntags_s):
            node_map3[nt] = ntags_d[i]

        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim != 3: break
            tag = vmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map3[i] = j

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map3[x] for x in item] for item in nodes]

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])            
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            done[3].append(tag)            
        return done, params

    # revolve
    # internally it is the same as extrude face
    def revolve_face_0D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve']=True
        return self.extrude_face_0D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)  
    def revolve_face_1D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve']=True
        return self.extrude_face_1D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)        
    def revolve_face_2D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve']=True
        return self.extrude_face_2D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)        
    def revolve_face_3D(self,  done, params, vdimtags, dimtags, dimtags2, *args, **kwargs):
        kwargs['revolve']=True
        return self.extrude_face_3D(done, params, vdimtags, dimtags, dimtags2, *args, **kwargs)        
        
    '''
    High-level interface test codes
    '''
    def extrude_surface_test(self, vtag, tag1, tag2, nlayers):
        
        from petram.geom.geom_utils import find_translate_between_surface
        from petram.geom.geom_utils import map_points_in_geom_info
        from petram.geom.geom_utils import map_lines_in_geom_info
        from petram.geom.geom_utils import map_surfaces_in_geom_info
        from petram.geom.geom_utils import map_volumes_in_geom_info        
        
        ptx, p, l, s, v = self.geom_info
        geom_data = (ptx, l, s, None)
        ax, an, px, d, affine, p_pairs, l_pairs = find_translate_between_surface(
                                                            [tag1], [tag2],
                                                            geom_data=geom_data)
        
        ents = gmsh.model.getBoundary(((2,tag1),), oriented=False)
        ents = list(ents)+[(2,tag1),]
        mdata = []
        for dim, tag in ents:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))
        
        ents = gmsh.model.getBoundary(((2,tag2),), oriented=False)        
        
        ws = self.prep_workspace()

        self.delete_all_except(2, [tag1])
        ret = gmsh.model.occ.extrude(((2,tag1),), d[0], d[1], d[2],
                                     numElements=[nlayers])
        gmsh.model.occ.synchronize()        

        info1 = self.geom_info
        info2 = self.read_geom_info()

        pmap, pmap_r = map_points_in_geom_info(info1, info2)
        lmap, lmap_r = map_lines_in_geom_info(info1, info2, pmap_r)
        smap, smap_r = map_surfaces_in_geom_info(info1, info2, lmap_r)
        vmap, vmap_r = map_volumes_in_geom_info(info1, info2, smap_r)
        
        self.show_all()
        gmsh.model.mesh.generate(1)
        ents = gmsh.model.getBoundary(((2,tag1),ret[0]), oriented=False)
        for dim, tag in ents:
            gmsh.model.mesh.setNodes(dim, tag, [], [], [])
            gmsh.model.mesh.setElements(dim, tag, [], [], [])

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        #print(node_map)
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # copy 1D elements on the source surface
        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim == 2: break
            tag = lmap[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
        print("node map here", node_map1)
        # copy 1D elements on the destination surface
        l_pairs2 = {x:lmap[l_pairs[x]] for x in l_pairs}
        R = affine[:3,:3]
        D = affine[:3,-1]
        node_map2 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info2[1][pmap[p_pairs[p]]]+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim == 2: break
            tag = l_pairs2[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map2[i] = j

            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]
            #print("setting", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        # copy source mesh
        for d in mdata:
            dim, tag, ndata, edata = d            
            if dim == 1: continue
            
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            
        # copy dest mesh
        for d in mdata:
            dim, tag, ndata, edata = d            
            if dim == 1: continue
            tag = ret[0][1]
            
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map2[i] = j

            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.hide(((2,tag1),ret[0]),)
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.generate(3)                    


        #
        #  copy back mash info
        #
        # first gather info
        ents = gmsh.model.getEntities()
        ents2 = gmsh.model.getBoundary(((2,tag1),), oriented=False)
        ents = [x for x in ents if not x in ents2 and x[0] != 0]
        ents = [x for x in ents if x[0] != 2 or x[1] != tag1]
        mdata = []
        for dim, tag in ents:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))

        
        # collect data for node mapping
        ents2_1D =  gmsh.model.getEntities(1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents2_1D]
        ntags_ws_1D = sum([x[0] for x in tmp], [])
        pos_ws_1D = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        ents2_2D = gmsh.model.getEntities(2)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents2_2D]
        ntags_ws_2D = sum([x[0] for x in tmp], [])
        pos_ws_2D = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        #node_map1 = {info2[1][k]+1:info1[1][pmap_r[k]]+1 for k in pmap_r}
        
        self.switch_model('main1')

        
        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}

        '''
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]
        print("nodemap1", node_map1)
        '''        
        #gmsh.model.occ.synchronize()
        #self.list_entities()                
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # send back 1D data
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 1: continue
            print("copy from ", dim, tag)                        
            tag = lmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            #print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}
        
        ents_1D = self.expand_dimtags(((3, vtag),), return_dim = 1)
        print("ents_1D", ents_1D)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]

        # send back 2D data
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 2: continue
            print("copy from ", dim, tag)            
            tag = smap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.show_all()            
        self.hide(((3, 2),), recursive=True)
        gmsh.model.mesh.generate(2)

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}

        ents_1D = self.expand_dimtags(((3, vtag),), return_dim = 1)                
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]

        ents_2D = self.expand_dimtags(((3, vtag),), return_dim = 2)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_2D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_2D]
        for i, nt in zip(idx, ntags_ws_2D):
            node_map1[nt] = ntags_m[i]
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 3: continue
            print("copy from ", dim, tag)            
            tag = vmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            #print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.show_all()            
        self.hide(((3, 2),), recursive=True)
        gmsh.model.mesh.generate(3)
            
    def run_generater(self, brep_input='', finalize=False, dim=3, msh_file='',
                      progressbar = None):
        
        if brep_input != '':
            filename = brep_input
        else:
            filename = self.target

        kwargs = {'CharacteristicLengthMax':self.clmax,
                  'CharacteristicLengthMin':self.clmin,
                  'EdgeResolution' : self.res,
                  'MeshAlgorithm'  : self.algorithm,
                  'MeshAlgorithm3D': self.algorithm3d,
                  'MaxThreads' : self.maxthreads}

        q = mp.Queue()
        p = mp.Process(target = generator,
                       args = (q, filename, self.mesh_sequence,
                               dim, finalize, msh_file, kwargs))
        p.start()
        istep = 0
        
        while True:
            try:
                ret = q.get(True, 1)
                if ret[0]: break
                if progressbar is not None:
                    istep += 1
                    progressbar.Update(istep, newmsg=ret[1])                    

            except QueueEmpty:
                if not p.is_alive():
                    assert False, "Child Process Died"
                    break
                time.sleep(1.)
                if progressbar is not None:
                    import wx
                    wx.Yield()
                    if progressbar.WasCancelled():
                       if p.is_alive():
                           p.terminate()
                       progressbar.Destroy()
                       assert False, "Mesh Generation Aborted"
                           
            time.sleep(0.01)
        return ret[1]

def generator(q, filename, sequence, dim, finalize, msh_file, kwargs):
    
    kwargs['queue'] = q
    mw = GMSHMeshWrapper(**kwargs)
    
    mw.mesh_sequence = sequence
    max_dim, done = mw.generate(dim=dim, brep_input = filename,
                                finalize=finalize, msh_file=msh_file)

    from petram.geom.read_gmsh import read_pts_groups, read_loops
        
    ptx, cells, cell_data = read_pts_groups(gmsh,
                                             finished_lines = done[1],
                                             finished_faces = done[2])
                
    data = ptx, cells, {}, cell_data, {}
    q.put((True, (max_dim, done, data)))
    
    

                   

        
