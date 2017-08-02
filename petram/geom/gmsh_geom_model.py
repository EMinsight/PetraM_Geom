'''

    Model Tree to stroe MFEM model parameters

'''
from __future__ import print_function

import tempfile
import os
import subprocess
import tempfile
import meshio
import numpy as np
import voropy

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshGeomModel')

from petram.model import Model

from petram.namespace_mixin import NS_mixin
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

debug = True

geom_key_dict = {'SurfaceBase': 'sb',
                 'Point': 'pt'}

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

class GmshPrimitiveBase(Model, Vtable_mixin, NS_mixin):
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
        return self.vt.panel_param(self)
        
    def get_panel1_value(self):
        return self.vt.get_panel_value(self)

    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return

    def import_panel1_value(self, v):
        return self.vt.import_panel_value(self, v)

    def panel1_tip(self):
        return self.vt.panel_tip()

    def build_geom(self, geom, objs):
        raise NotImplementedError(
             "you must specify this method in subclass")
    
    def get_special_menu(self):
        return [('Build this step', self.onBuildAfter)]

    def onBuildAfter(self, evt):
        dlg = evt.GetEventObject()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        
        self.parent.build_geom(stop2 = self)
        dlg.OnRefreshTree()        
        evt.Skip()


class GmshGeom(Model, NS_mixin):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(GmshGeom, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    def get_possible_child(self):
        from .gmsh_primitives import Circle, Rect, Polygon, Extrude, Revolve, Difference
        return [Circle, Rect, Polygon, Extrude, Revolve, Difference]
    
    def get_special_menu(self):
        return [('Build all', self.onBuildAll)]

    def onBuildAll(self, evt):
        dlg = evt.GetEventObject()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        
        self.build_geom()
        dlg.OnRefreshTree()
        evt.Skip()
        
    def build_geom(self, stop1=None, stop2=None):
        children = [x for x in self.walk()]
        children = children[1:]

        objs = GeomObjs()
        
        import pygmsh
        geom = pygmsh.Geometry()
        for child in children:
            if not child.enabled: continue            
            child.vt.preprocess_params(child)
            if child is stop1: break            # for build before
            child.build_geom(geom, objs)
            if child is stop2: break            # for build after

        txt_unrolled, num_entities = generate_mesh(geom, dim = 0)
        if debug:
            for l in txt_unrolled:
                 print(l.strip())
        self._txt_unrolled = [x.strip() for x in txt_unrolled]
        self._num_entities = num_entities
        print(num_entities)
        self._objs = objs

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
            fid.write('\n'.join(self._txt_unrolled))
            fid.close()
        
    def panel1_param(self):
        return [["", "Geometry model using GMSH", 2, None],
                [None, None, 141, {"label": "Export...",
                                   "func": self.onExportGeom,
                                   "noexpand": True}],]
                
    def get_panel1_value(self):
        return [None, None]
       
    def import_panel1_value(self, v):
        pass

def generate_mesh(
        geo_object,
        optimize=True,
        num_quad_lloyd_steps=10,
        num_lloyd_steps=1000,
        verbose=True,
        dim=3,
        prune_vertices=True
        ):
    from pygmsh.helper import _get_gmsh_exe
    
    handle, geo_filename = tempfile.mkstemp(suffix='.geo')
    os.write(handle, geo_object.get_code().encode())
    if dim == 0:
        os.write(handle, "Geometry.OldNewReg=0;\n")
        os.write(handle, 'Printf("Number of entitites, : %g, %g, %g, %g :", newp, newl, news, newv);\n')
    os.close(handle)

    gmsh_executable = _get_gmsh_exe()

    if dim > 0:
        handle, msh_filename = tempfile.mkstemp(suffix='.msh')
        os.close(handle)
        cmd = [
            gmsh_executable,
            '-{}'.format(dim), '-bin', geo_filename, '-o', msh_filename
            ]
        gmsh_major_version = geo_object.get_gmsh_major()
        if gmsh_major_version < 3 and optimize:
            cmd += ['-optimize']
        if num_quad_lloyd_steps > 0:
            cmd += ['-optimize_lloyd', str(num_quad_lloyd_steps)]
        
    else:
        handle, geou_filename = tempfile.mkstemp(suffix='.geo_unrolled')
        os.close(handle)
        cmd = [
            gmsh_executable,
            '-{}'.format(dim), geo_filename, '-o', geou_filename
            ]

    # http://stackoverflow.com/a/803421/353337
    p = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
    stdoutdata = p.stdout.readlines()    
    if verbose:
        for x in stdoutdata:
            print(x.decode('utf-8'), end='')

    p.communicate()
    assert p.returncode == 0,\
        'Gmsh exited with error (return code {}).'.format(p.returncode)
    
    if dim == 0:
        for x in stdoutdata:
            if x.find("Number of entitites,")!= -1:
                txt = x[x.find("Number of entitites,"):].split(":")
                print(txt)
                num_entities = [int(x)-1 for x in txt[1].split(',')]
                
        fid = open(geou_filename, 'r')
        lines = fid.readlines()
        fid.close()
        if not debug:
            os.remove(geo_filename)
            os.remove(geou_filename)
        return lines, num_entities

    X, cells, pt_data, cell_data, field_data = meshio.read(msh_filename)

    # clean up
    if not debug:    
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
    
