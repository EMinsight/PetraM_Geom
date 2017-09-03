'''

    Model Tree to stroe MFEM model parameters

'''
from __future__ import print_function

import tempfile
import os
import subprocess
import traceback
import sys

import numpy as np
import warnings

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshGeomModel')

from petram.model import Model

from petram.geom.geom_model import GeomBase
from petram.namespace_mixin import NS_mixin
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

debug = True

geom_key_dict = {'SurfaceBase': 'sb',
                 'PlaneSurface' : 'sp',
                 'Point': 'pt',
                 'Line': 'ln',
                 'Spline': 'sp'}

def get_geom_key(obj):
    if obj.__class__ in geom_key_dict:
        return geom_key_dict[obj.__class__]
    name = obj.__class__.__name__
    key = ''.join([i.lower() for i in name if not i.isupper()])

    for k in geom_key_dict.keys():
        if geom_key_dict[k] == key:
            assert False, key + " is used for " + k.__name__
            
    geom_key_dict[obj.__class__] = key
    if debug: print(geom_key_dict)
    return key

class GeomObjs(dict):
    def addobj(self, obj, name):
        key = ''.join([i for i in name if not i.isdigit()])        
        keys = self.keys()
        nums = []
        for k in keys:
           t = ''.join([i for i in k if not i.isdigit()])
           if t == key:
              n = int(''.join([i for i in k if i.isdigit()]))
              nums.append(n)
        if len(nums) == 0:
           newkey = key+str(1)
        else:
           newkey = key+str(max(nums)+1)
        self[newkey] = obj
        return newkey

class GmshPrimitiveBase(GeomBase, Vtable_mixin):
    hide_ns_menu = True
    has_2nd_panel = False
    isGeom = True
        
    def __init__(self, *args, **kwargs):
        super(GmshPrimitiveBase, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    def attribute_set(self, v):
        v = super(GmshPrimitiveBase, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v
    
    def panel1_param(self):
        from wx import BU_EXACTFIT
        b1 = {"label": "S", "func": self.onBuildBefore,
              "noexpand": True, "style": BU_EXACTFIT}
        b2 = {"label": "R", "func": self.onBuildAfter,
              "noexpand": True, "style": BU_EXACTFIT}
        
        ll = [[None, None, 241, {'buttons':[b1,b2],
                                 'alignright':True,
                                 'noexpand': True},],]
        ll.extend(self.vt.panel_param(self))
        return ll
        
    def get_panel1_value(self):
        return [None] + list(self.vt.get_panel_value(self))

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        return self.vt.import_panel_value(self, v[1:])

    def panel1_tip(self):
        return [None]+list(self.vt.panel_tip())

    def build_geom(self, geom, objs):
        raise NotImplementedError(
             "you must specify this method in subclass")

    def gsize_hint(self, geom, objs):
        '''
        return quick estimate of geometry size min and max
        '''
        warnings.warn("Not implemented", Warning)
        return -1, -1
    
    def get_special_menu(self):
        return [('Build this step', self.onBuildAfter)]

    def _onBuildThis(self, evt, **kwargs):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        try:
            self.parent.build_geom(**kwargs)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to build geometry',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.parent.onUpdateGeoView(evt)

        
    def onBuildBefore(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        mm = dlg.get_selected_mm()
        self._onBuildThis(evt, stop1 = mm)
        evt.Skip()
        
    def onBuildAfter(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        mm = dlg.get_selected_mm()
        self._onBuildThis(evt, stop2 = mm)
        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.select_next_enabled()
        evt.Skip()
        
class GmshGeom(GeomBase):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(GmshGeom, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    def get_possible_child(self):
        from .gmsh_primitives import Point, Line, Spline, Circle, Rect, Polygon, Extrude, Revolve, LineLoop, CreateLine, CreateSurface, CreateVolume, SurfaceLoop, Union, Intersection, Difference, Fragments
        return [Point, Line, Circle, Rect, Polygon, Spline, CreateLine, CreateSurface, CreateVolume, LineLoop, SurfaceLoop, Extrude, Revolve, Union, Intersection, Difference, Fragments]
    
    def get_special_menu(self):
        return [('Build All', self.onBuildAll),
                ('Export .geo', self.onExportGeom)]
    
    def panel1_param(self):
        return [["", "Geometry model using GMSH", 2, None],
                [None, None, 341, {"label": "Build All",
                                   "func": 'onBuildAll',
                                   "noexpand": True}],]
                
    def get_panel1_value(self):
        return [None, self]
       
    def import_panel1_value(self, v):
        pass
    

    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()

        try:
            self.build_geom(finalize = True)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to build geometry',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()

        filename = os.path.join(viewer.model.owndir(), self.name())
        self.onUpdateGeoView(evt, filename = filename)
        fid = open(filename + '.geo_unrolled', 'w')
        fid.write('\n'.join(self._txt_unrolled))
        fid.close()
        evt.Skip()
        
    def onUpdateGeoView(self, evt, filename = None):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        geo_text = self._txt_rolled[:]
        xyz = guess_geom_size(self._txt_unrolled)

        clmax = np.max(np.max(xyz, 0) - np.min(xyz, 0))/3.
        clmin = np.min(np.max(xyz, 0) - np.min(xyz, 0))/3.        
        self._clmax_guess = (clmax, clmin)
        geo_text.extend(['Show "*";',
                         'Mesh.CharacteristicLengthMax = '+str(clmax) + ';'])

        ret =  generate_mesh(geo_object = None,
                             dim = 2,
                             filename = filename,
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text)
        viewer.set_figure_data('geom', self.name(), ret)
        viewer.update_figure('geom', self.name())

        self._geom_coords = ret
        viewer._s_v_loop['geom'] = read_loops(self._txt_unrolled)
        viewer._s_v_loop['mesh'] = viewer._s_v_loop['geom']
        
    def build_geom(self, stop1=None, stop2=None, filename = None,
                   finalize = False):
        '''
        filename : export geometry to a real file (for debug)
        '''
        children = [x for x in self.walk()]
        children = children[1:]

        objs = GeomObjs()
        self._objs = objs        
        import pygmsh
        geom = pygmsh.Geometry()
        geom.set_factory('OpenCASCADE')
        
        for child in children:
            if not child.enabled: continue            
            child.vt.preprocess_params(child)
            if child is stop1: break            # for build before
            child.build_geom(geom, objs)
            if child is stop2: break            # for build after

        unrolled, rolled, entities = generate_mesh(geom, dim = 0,
                                                  filename = filename)
        unrolled = [x.strip() for x in unrolled]
        s, v =  read_loops(unrolled)
        if finalize:
            dim = check_dim(unrolled)
            extra = []
            if ((dim == 2 and len(s.keys()) > 1) or
                (dim == 3 and len(s.keys()) > 1 and len(v.keys()) == 1)):
                print("splitting surface",  s.keys())
                extra.append(BoolFramgents_extra('final_s', 'Surface', s.keys()))
                unrolled, rolled, entities = generate_mesh(geom, dim = 0,
                                                             filename = filename,
                                                             extra = extra)
                unrolled = [x.strip() for x in unrolled]
                s, v =  read_loops(unrolled)              
            if dim == 3 and len(v.keys()) > 1:                
                print("splitting volume",  v.keys())                
                extra.append(BoolFramgents_extra('final_v', 'Volume', v.keys()))
                unrolled, rolled, entities = generate_mesh(geom, dim = 0,
                                                           filename = filename,
                                                           extra = extra)
                unrolled = [x.strip() for x in unrolled]
                s, v =  read_loops(unrolled)
        self._txt_rolled = rolled                
        self._txt_unrolled = unrolled
        self._num_entities = entities


    def onExportGeom(self, evt):
        print("export geom file")
        if not hasattr(self, "_txt_unrolled"):
            evt.Skip()
            return
        from ifigure.widgets.dialog import write
        parent = evt.GetEventObject()
        path = write(parent,
                     message = 'Enter .geo file name',
                     wildcard = '*.geo')
        if path != '':
            fid = open(path, 'w')
            fid.write('\n'.join(self._txt_rolled))
            fid.close()
         
    def load_gui_figure_data(self, viewer):
        import meshio
        filename = os.path.join(viewer.model.owndir(), self.name())
        msh_filename = filename + '.msh'
        ret = meshio.read(msh_filename)

        filename = os.path.join(viewer.model.owndir(), self.name())
        filename = filename + '.geo_unrolled'
        if os.path.exists(filename):
            fid = open(filename, 'r')
            unrolled = [l.strip() for l in fid.readlines()]
            fid.close()
            viewer._s_v_loop['geom'] = read_loops(unrolled)        
            viewer._s_v_loop['mesh'] = viewer._s_v_loop['geom']

        return 'geom', self.name(), ret
    
    def is_viewmode_grouphead(self):
        return True
    
   
def check_dim(unrolled):
    for line in unrolled:
        if line.startswith('Volume'): return 3
    for line in unrolled:
        if line.find('Surface'): return 2
    return 1

def guess_geom_size(unrolled):
    points = []
    for line in unrolled:
        if line.startswith("Point("):
            try:
                coords = line.split("=")[1]
                coords = coords[coords.find("{")+1:coords.find("}")]
                xyz = np.array([float(x) for x in coords.split(",")[:3]])
                points.append(xyz)
            except:
                pass
    points = np.vstack(points)
    return points
            
def read_loops(unrolled):
    ll = {}  # line loop
    sl = {}  # surface loop
    v = {}
    s = {}

    def split_line(line):
        a, b = line.split("=")
        k = int(a[a.find("(")+1:a.find(")")])
        idx = [abs(int(x)) for x in b[b.find("{")+1:b.find("}")].split(",")]
        return k, idx
    
    for line in unrolled:
        line = line.strip()
        if line.startswith("Surface Loop("):
            k, idx = split_line(line)
            sl[k] = idx
        elif line.startswith("Line Loop("):
            k, idx = split_line(line)
            ll[k] = idx
        elif line.startswith("Volume("):
            k, idx = split_line(line)
            v[k] = idx
        elif line.startswith("Plane Surface("):
            k, idx = split_line(line)
            s[k] = idx
        elif line.startswith("Surface("):
            k, idx = split_line(line)
            s[k] = idx
        else:
            pass

    for kv in v.keys():
        tmp = []
        for k in v[kv]:
            tmp.extend(sl[k])
        v[kv] = list(set(tmp))
    for ks in s.keys():
        tmp = []
        for k in s[ks]:
            tmp.extend(ll[k])
        s[ks] = list(set(tmp))
    return s, v

def BoolFramgents_extra(name, shape_type,  inputs, delete = True):
    txt = '{}[] = {}{{{} {{{}}}; {}}} {{{} {{{}}}; {}}};'.format(
                name, 
                'BooleanFragments',
                shape_type,
                ','.join(str(e) for e in inputs[:1]),
                'Delete;' if delete else '',
                shape_type,
                ','.join(str(e) for e in inputs[1:]),
                'Delete;' if delete else '')
    return txt

def generate_mesh(
        geo_object = None,
        optimize=True,
        num_quad_lloyd_steps=10,
        num_lloyd_steps=1000,
        verbose=True,
        dim=3,
        prune_vertices=True,
        filename=None,
        geo_text=None,
        extra=None,
        bin = '-bin',
        verbosity = '4'
        ):
    from pygmsh.helper import _get_gmsh_exe, _is_flat
    import meshio
    import voropy    

    if filename is None:
        handle, geo_filename = tempfile.mkstemp(suffix='.geo')
    else:
        geo_filename = filename + '.geo'
        handle = os.open(geo_filename,  os.O_WRONLY | os.O_CREAT |
                         os.O_TRUNC)
    extra = [] if extra is None else extra    

    if geo_object is not None:
       rolled =  geo_object.get_code().encode().split("\n")
    elif geo_text is not None:
       rolled = geo_text

    if dim == 0:
        rolled.extend(extra)
        rolled.append("Geometry.OldNewReg=0;\n")
        rolled.append('Printf("Number of entitites, : %g, %g, %g, %g :", newp, newl, news, newv);\n')
    os.write(handle, "\n".join(rolled))
    os.close(handle)

    gmsh_executable = _get_gmsh_exe()

    if dim > 0:
        if filename is None:
            handle, msh_filename = tempfile.mkstemp(suffix='.msh')
            os.close(handle)            
        else:
            msh_filename = filename + '.msh'
        cmd = [
            gmsh_executable, '-v', verbosity, 
            '-{}'.format(dim), bin, geo_filename, '-o', msh_filename
            ]
        if num_quad_lloyd_steps > 0:
            cmd += ['-optimize_lloyd', str(num_quad_lloyd_steps)]
    elif dim < 0:
        if filename is None:
            handle, msh_filename = tempfile.mkstemp(suffix='.msh')
            os.close(handle)            
        else:
            msh_filename = filename + '.msh'
        cmd = [
            gmsh_executable, bin, geo_filename, '-o', msh_filename
            ]
    else:
        if filename is None:
            handle, geou_filename = tempfile.mkstemp(suffix='.geo_unrolled')
            os.close(handle)         
        else:
            geou_filename = filename + '.geo_unrolled'
        cmd = [
            gmsh_executable,
            '-{}'.format(dim), geo_filename, '-o', geou_filename
            ]
    cmd = [x for x in cmd if x != '']
    # http://stackoverflow.com/a/803421/353337
    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        bufsize = 0)
    
    stdoutdata = []
    for line in iter(p.stdout.readline, ''):
        if verbose:
            print(line.decode('utf-8'), end='')
            sys.stdout.flush()
        stdoutdata.append(line)

    p.communicate()
    assert p.returncode == 0,\
        'Gmsh exited with error (return code {}).'.format(p.returncode)
    
    if dim == 0:
        for x in stdoutdata:
            if x.find("Number of entitites,")!= -1:
                txt = x[x.find("Number of entitites,"):].split(":")
                num_entities = [int(x)-1 for x in txt[1].split(',')]
                
        fid = open(geou_filename, 'r')
        lines = fid.readlines()
        fid.close()
        #if filename is None:            
        #    os.remove(geo_filename)
        #    os.remove(geou_filename)
        return lines, rolled, num_entities

    # meshio does not read $Periodic....
    fid = open(msh_filename, 'r')
    lines = fid.readlines()
    fid.close()
    has_periodic = False
    for n1, l in enumerate(lines):
        if l.strip() == '$Periodic':
            has_periodic = True
            ps = n1
            break
    if has_periodic:
        for n2, l in enumerate(lines):
            if l.strip() == '$EndPeriodic':
                pe = n2
                break
        print(n1, n2)
        lines = lines[:ps]+lines[pe+1:]
        fid = open(msh_filename, 'w')
        fid.write(''.join(lines))
        fid.close()
        
    X, cells, pt_data, cell_data, field_data = meshio.read(msh_filename)

    # clean up
    if filename is None:    
        os.remove(geo_filename)
        os.remove(msh_filename)

    # Lloyd smoothing
    if not _is_flat(X) or 'triangle' not in cells:
        if verbose:
            print(
                'Not performing Lloyd smoothing '
                '(only works for flat triangular meshes).'
                )
        return X, cells, pt_data, cell_data, field_data
    if dim < 0:
        return X, cells, pt_data, cell_data, field_data
    
    if num_lloyd_steps == 0 and num_quad_lloyd_steps == 0:
        return X, cells, pt_data, cell_data, field_data
    
    if verbose:
        print('Lloyd smoothing...')
    # find submeshes
    a = cell_data['triangle']['geometrical']
    # http://stackoverflow.com/q/42740483/353337
    submesh_bools = {v: v == a for v in np.unique(a)}

    X, cells['triangle'] = voropy.smoothing.lloyd_submesh(
            X, cells['triangle'], submesh_bools,
            tol=0.0, max_steps=num_lloyd_steps,
            verbose=False
            )

    if prune_vertices:
        # Make sure to include only those vertices which belong to a triangle.
        uvertices, uidx = np.unique(cells['triangle'], return_inverse=True)
        cells = {'triangle': uidx.reshape(cells['triangle'].shape)}
        cell_data = {'triangle': cell_data['triangle']}
        X = X[uvertices]
        for key in pt_data:
            pt_data[key] = pt_data[key][uvertices]

    return X, cells, pt_data, cell_data, field_data
    
