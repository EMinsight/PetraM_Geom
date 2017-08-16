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
import traceback
import warnings

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
        self._onBuildThis(evt, stop1 = self)
        evt.Skip()
        
    def onBuildAfter(self, evt):        
        self._onBuildThis(evt, stop2 = self)
        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.select_next_enabled()
        evt.Skip()
        
    def onItemSelChanged(self, evt):
        '''
        GUI response when model object is selected in
        the dlg_edit_model
        '''
        viewer = evt.GetEventObject().GetTopLevelParent().GetParent()
        viewer.canvas.use_navibar_palette('petram_geom', mode = '3D')
                
class GmshGeom(Model, NS_mixin):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(GmshGeom, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    def get_possible_child(self):
        from .gmsh_primitives import Circle, Rect, Polygon, Extrude, Revolve, Difference
        return [Circle, Rect, Polygon, Extrude, Revolve, Difference]
    
    def get_special_menu(self):
        return [('Build All', self.onBuildAll)]

    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()

        try:
            self.build_geom()
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to build geometry',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.onUpdateGeoView(evt)
        evt.Skip()
        
    def onUpdateGeoView(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        geo_text = self._txt_unrolled[:]
#        geo_text.extend(['Show "*";',
#                         'Transfinite Line "*"  = 10;'])
        geo_text.extend(['Show "*";',
                         'Mesh.CharacteristicLengthMax = 0.1;'])

        ret =  generate_mesh(geo_object = None,
                             dim = 2,
#                             filename = "/Users/shiraiwa/test",
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text)
        from .geo_plot import plot_geometry
        plot_geometry(viewer, ret)
    
    def build_geom(self, stop1=None, stop2=None):
        children = [x for x in self.walk()]
        children = children[1:]

        objs = GeomObjs()
        self._objs = objs        
        import pygmsh
        geom = pygmsh.Geometry()
        for child in children:
            if not child.enabled: continue            
            child.vt.preprocess_params(child)
            if child is stop1: break            # for build before
            child.build_geom(geom, objs)
            if child is stop2: break            # for build after

        txt_unrolled, num_entities = generate_mesh(geom, dim = 0)

        self._txt_unrolled = [x.strip() for x in txt_unrolled]
        self._num_entities = num_entities


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
                [None, None, 141, {"label": "Build All",
                                   "func": self.onBuildAll,
                                   "noexpand": True}],
                [None, None, 141, {"label": "Export...",
                                   "func": self.onExportGeom,
                                   "noexpand": True}],]
                
    def get_panel1_value(self):
        return [None, None, None]
       
    def import_panel1_value(self, v):
        pass
    
    def onItemSelChanged(self, evt):
        '''
        GUI response when model object is selected in
        the dlg_edit_model
        '''
        viewer = evt.GetEventObject().GetTopLevelParent().GetParent()
        viewer.canvas.use_navibar_palette('petram_geom', mode = '3D')
        
def generate_mesh(
        geo_object = None,
        optimize=True,
        num_quad_lloyd_steps=10,
        num_lloyd_steps=1000,
        verbose=True,
        dim=3,
        prune_vertices=True,
        filename = None,
        geo_text = None,
        ):
    from pygmsh.helper import _get_gmsh_exe, _is_flat

    if filename is None:
        handle, geo_filename = tempfile.mkstemp(suffix='.geo')
    else:
        geo_filename = filename + '.geo'
        handle = os.open(geo_filename,  os.O_WRONLY | os.O_CREAT |
                         os.O_TRUNC)
        

    if geo_object is not None:
       os.write(handle, geo_object.get_code().encode())
    elif geo_text is not None:
       os.write(handle, '\n'.join(geo_text))

    if dim == 0:
        os.write(handle, "Geometry.OldNewReg=0;\n")
        os.write(handle, 'Printf("Number of entitites, : %g, %g, %g, %g :", newp, newl, news, newv);\n')
    os.close(handle)

    gmsh_executable = _get_gmsh_exe()

    if dim > 0:
        if filename is None:
            handle, msh_filename = tempfile.mkstemp(suffix='.msh')
            os.close(handle)            
        else:
            msh_filename = filename + '.msh'
        cmd = [
            gmsh_executable,
            '-{}'.format(dim), '-bin', geo_filename, '-o', msh_filename
            ]
        if num_quad_lloyd_steps > 0:
            cmd += ['-optimize_lloyd', str(num_quad_lloyd_steps)]
        
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
                num_entities = [int(x)-1 for x in txt[1].split(',')]
                
        fid = open(geou_filename, 'r')
        lines = fid.readlines()
        fid.close()
        if filename is None:            
            os.remove(geo_filename)
            os.remove(geou_filename)
        return lines, num_entities

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
    
