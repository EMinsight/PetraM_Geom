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

debug = False
        
class GmshMeshActionBase(Mesh, Vtable_mixin):
    hide_ns_menu = True
    has_2nd_panel = False
    isGmshMesh = True
    
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
        v['algorithm'] = 'default'
        v['algorithm3d'] = 'default'                
        super(GmshMesh, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v
    
    def panel1_param(self):    
        ll =   [["Geometry", self.geom_group,  0, {},],]
        ll.extend(self.vt.panel_param(self))

        from .gmsh_mesher import MeshAlgorithm, MeshAlgorithm3D

        c1 = MeshAlgorithm.keys()
        c2 = MeshAlgorithm3D.keys()

        from wx import CB_READONLY
        setting1 = {"style":CB_READONLY, "choices": c1}
        setting2  ={"style":CB_READONLY, "choices": c2}
        
        ll.append(["2D Algorithm", c1[-1], 4, setting1])
        ll.append(["3D Algorithm", c2[-1], 4, setting2])
        ll.append([None, None, 141, {"label": "Use default size",
                                  "func": self.onSetDefSize,
                                   "noexpand": True}])
        ll.append([None, None, 141, {"label": "Build All",
                                  "func": self.onBuildAll,
                                   "noexpand": True}])
        return ll
    
    def get_panel1_value(self):
        return ([self.geom_group,] + list(self.vt.get_panel_value(self)) +
                [self.algorithm, self.algorithm3d, None, None])
    
    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return
    
    def import_panel1_value(self, v):
        self.geom_group = str(v[0])        
        self.vt.import_panel_value(self, v[1:-4])

        from .gmsh_mesher import MeshAlgorithm, MeshAlgorithm3D

        self.algorithm = str(v[-4])
        self.algorithm3d = str(v[-3])
        
    def panel1_tip(self):
        return [None] + self.vt.panel_tip() + [None]*4
        
    def get_possible_child(self):
        from .gmsh_mesh_actions import TransfiniteLine, FreeFace, FreeVolume, FreeEdge, CharacteristicLength, Rotate, Translate, CopyFace
        return [FreeVolume, FreeFace, FreeEdge, TransfiniteLine, CharacteristicLength, Rotate, Translate, CopyFace]
    
    def get_special_menu(self):
        return [('Build All', self.onBuildAll),
                ('Export .geo', self.onExportGeom),
                ('Export .msh', self.onExportMsh)]    
    
    def onSetDefSize(self, evt):
        geom_root = self.root()['Geometry'][self.geom_group]                
        clmax_root, clmin_root = geom_root._clmax_guess
        self.clmax_txt = str(clmax_root)
        self.clmin_txt = str(clmin_root)        
        dlg = evt.GetEventObject().GetTopLevelParent()        
        dlg.OnItemSelChanged()
        
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
        geo_text = self.assign_physical(self._txt_rolled[:])
        if path != '':
            fid = open(path, 'w')
            fid.write('\n'.join(geo_text))
            fid.close()
            
    def onExportMsh(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        src = os.path.join(viewer.model.owndir(), self.name())+'.msh'

        from ifigure.widgets.dialog import write
        parent = evt.GetEventObject()        
        dst = write(parent,
                     message = 'Enter .msh file name',
                     wildcard = '*.msh')
        if dst == '': return
        try:
            import shutil
            shutil.copyfile(src, dst)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to export msh file',
                                 title='Error',
                                 traceback=traceback.format_exc())
        
            
    def onUpdateMeshView(self, evt, filename=None, bin='-bin', geo_text = None):
        from petram.geom.gmsh_geom_model import read_loops, generate_mesh
        from petram.geom.geo_plot import plot_geometry, oplot_meshed        
        
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        geo_text = self._txt_rolled[:] if geo_text is None else geo_text
        
        ret =  generate_mesh(geo_object = None,
                             dim = self._max_mdim,
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text,
                             filename=filename, bin=bin)

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
        
        filename = os.path.join(viewer.model.owndir(), self.name())
        geo_text = self.assign_physical(self._txt_rolled[:])        
        self.onUpdateMeshView(evt, filename = filename, bin='',
                              geo_text = geo_text)

        evt.Skip()

    def assign_physical(self, geo_text):
        from petram.geom.gmsh_geom_model import read_loops
        geom_root = self.root()['Geometry'][self.geom_group]        
        s, v = read_loops(geom_root._txt_unrolled)
        

        has_finalv = False
        has_finals = False
        has_finall = False
        for line in geo_text:
            if line.startswith('final_v[]'): has_finalv = True
            if line.startswith('final_s[]'): has_finals = True
            if line.startswith('final_l[]'): has_finall = True

        t1 = 'final_s() = Unique(Abs(Boundary{ Volume{final_v()}; }));'
        t2 = 'final_l() = Unique(Abs(Boundary{ Surface{final_s()}; }));'
        t3 = 'final_p() = Unique(Abs(Boundary{ Line{final_l()}; }));'
        
        tt1= ['ipv = 0;',
              'For ii In {0 : #final_v[]-1}',
              '   ipv = ipv+1;',
              '   Physical Volume (StrCat("volume", Sprintf("%g", final_v[ii])),ipv) = {final_v[ii]};',
              'EndFor',]
        tt2 = ['ips = 0;',
               'For ii In {0 : #final_s[]-1}',
               '   ips = ips+1;',
               '   Physical Surface (StrCat("surface", Sprintf("%g", final_s[ii])),ips) = {final_s[ii]};',
               'EndFor']
        tt3 = ['ipl = 0;',
               'For ii In {0 : #final_l[]-1}',
               '   ipl = ipl+1;',
               '   Physical Line (StrCat("line", Sprintf("%g", final_l[ii])),ipl) = {final_l[ii]};',
              'EndFor',]

        ## if volume (surface loop) exitsts but there is only one volume,
        ## set it to final_v
        if len(v.keys()) == 1:
            ipv = str(v.keys()[0])
            geo_text.append('final_v[] = {'+ipv+'};')
            has_finalv = True
            has_finals = False
        ## if surface (line loop)  exitsts but there is only one volume,
        ## set it to final_s
        if len(s.keys()) == 1:
            ips = str(s.keys()[0])
            geo_text.append('final_s[] = {'+ips+'};')
            has_finals = True
        if has_finalv:
            geo_text.extend([t1, t2, t3]+ tt1 + tt2 + tt3)
        if has_finals:
            geo_text.extend([t2, t3] + tt2 + tt3)
        if has_finall:
            geo_text.extend([t3] + tt3)
        return geo_text
            
            
    def build_mesh(self, geom_root, stop1=None, stop2=None, filename = None):
        self.vt.preprocess_params(self)
        
        num_entities = geom_root._num_entities
        geom_coords = geom_root._geom_coords
        children = [x for x in self.walk()]
        children = children[1:]

        from .gmsh_mesher import GmshMesher
        mesher = GmshMesher(num_entities,
                            geom_coords = geom_coords,
                            CharacteristicLengthMax = self.clmax,
                            CharacteristicLengthMin = self.clmin,
                            MeshAlgorithm = self.algorithm,
                            MeshAlgorithm3D = self.algorithm3d)                


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
        

        
