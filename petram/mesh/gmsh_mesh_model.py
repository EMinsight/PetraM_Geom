from __future__ import print_function

import tempfile
import os
import subprocess
import tempfile
import weakref
import numpy as np
import traceback

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomModel')

from petram.mesh.mesh_model import Mesh
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

debug = True

        
class GmshMeshActionBase(Mesh, Vtable_mixin):
    hide_ns_menu = True
    has_2nd_panel = False
    isGmshMesh = True
    
#    def __init__(self, *args, **kwargs):
#        super(GmshMeshActionBase, self).__init__(*args, **kwargs)
#        NS_mixin.__init__(self, *args, **kwargs)

    def attribute_set(self, v):
        v = super(GmshMeshActionBase, self).attribute_set(v)
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
        self.vt.import_panel_value(self, v[1:])
        return True

    def panel1_tip(self):
        return [None] + self.vt.panel_tip()

    def add_meshcommand(self):
        raise NotImplementedError(
             "you must specify this method in subclass")
    
    def _onBuildThis(self, evt, **kwargs):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()        
        geom_root = self.root()['Geometry'][self.parent.geom_group]
        try:
            self.parent.build_mesh(geom_root, **kwargs)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate meshing script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.parent.onUpdateMeshView(evt)
        
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

    def element_selection_empty(self):
        return {'volume':[],
                'face':[],
                'edge':[],
                'point':[],}, None
    
    def get_element_selection(self):
        return self.element_selection_empty()
                
    def onItemSelChanged(self, evt):
        super(GmshMeshActionBase, self).onItemSelChanged(evt)
        dlg = evt.GetEventObject().GetTopLevelParent()
        self.update_viewer_selection(dlg)
        
    def update_after_ELChanged(self, dlg):
        self.update_viewer_selection(dlg)
        
    def update_viewer_selection(self, dlg):
        viewer = dlg.GetParent()                
        sel, mode = self.get_element_selection()
        figobj = viewer.highlight_element(sel)

        viewer.set_sel_mode(mode)
        viewer.set_sel_mode() # update buttons        
        if figobj is not None:
            import ifigure.events
            sel = [weakref.ref(figobj._artists[0])]
            ifigure.events.SendSelectionEvent(figobj, dlg, sel)                
            

data = (('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size(def)',
                                default = 1.0, 
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size(def)',
                                default = 1.0, 
                                tip = "CharacteristicLengthMin" )),)
                
class GmshMesh(Mesh, Vtable_mixin):
    has_2nd_panel = False
    isMeshGroup = True
    vt = Vtable(data)    
#    def __init__(self, *args, **kwargs):
#        super(GmshMesh, self).__init__(*args, **kwargs)
#        NS_mixin.__init__(self, *args, **kwargs)
        
    def attribute_set(self, v):
        v['geom_group'] = ''
        super(GmshMesh, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v
    
    def panel1_param(self):    
        ll =   [["Geometry", self.geom_group,  0, {},],]
        ll.extend(self.vt.panel_param(self))        
        ll.append([None, None, 141, {"label": "Build All",
                                  "func": self.onBuildAll,
                                   "noexpand": True}])

        return ll
    
    def get_panel1_value(self):
        return [self.geom_group,] + list(self.vt.get_panel_value(self)) + [None]
    
    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return
    
    def import_panel1_value(self, v):
        self.geom_group = str(v[0])        
        self.vt.import_panel_value(self, v[1:])
        
    def panel1_tip(self):
        return [None] + self.vt.panel_tip() + [None]
        
    def get_possible_child(self):
        from .gmsh_mesh_actions import TransfiniteLine, FreeFace, FreeVolume, FreeEdge, CharacteristicLength, Rotate, Translate
        return [FreeVolume, FreeFace, FreeEdge, TransfiniteLine, CharacteristicLength, Rotate, Translate]
    
    def get_special_menu(self):
        return [('Build All', self.onBuildAll),
                ('Export .geo', self.onExportGeom)]

    def onExportGeom(self, evt):
        print("export geom file")
        if not hasattr(self, "_txt_rolled"):
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
            
    def onUpdateMeshView(self, evt):
        from petram.geom.gmsh_geom_model import read_loops, generate_mesh
        from petram.geom.geo_plot import plot_geometry, oplot_meshed        
        
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        geo_text = self._txt_rolled[:]        
        ret =  generate_mesh(geo_object = None,
                             dim = self._max_mdim,
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text)

        oplot_meshed(viewer, ret)        
        geom_root = self.root()['Geometry'][self.geom_group]        
        viewer._s_v_loop = read_loops(geom_root._txt_unrolled)
    
    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        
        geom_root = self.root()['Geometry'][self.geom_group]
        if not hasattr(geom_root, "_txt_unrolled"):
            geom_root.onBuildAll(evt)
            
        try:
           self.build_mesh(geom_root)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate mesh script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.onUpdateMeshView(evt)        
        evt.Skip()
        
    def build_mesh(self, geom_root, stop1=None, stop2=None):
        lines = [x.strip() for x in geom_root._txt_unrolled]
        num_entities = geom_root._num_entities
        geom_coords = geom_root._geom_coords
        children = [x for x in self.walk()]
        children = children[1:]

        clmax_root = geom_root._clmax_guess
        clmin_root = clmax_root/100.
        from .gmsh_mesher import GmshMesher
        mesher = GmshMesher(num_entities,
                            geom_coords = geom_coords,
                            CharacteristicLengthMax = clmax_root,
                            CharacteristicLengthMin = clmin_root)


        for child in children:
            if not child.enabled: continue            
            child.vt.preprocess_params(child)
            if child is stop1: break            # for build before
            child.add_meshcommand(mesher)
            if child is stop2: break            # for build after
            
        lines, max_mdim = mesher.generate()
        if debug:
            for l in lines:
                 print(l)

        lines = geom_root._txt_rolled + lines
        self._txt_rolled = lines
        self._max_mdim = max_mdim
        

        
