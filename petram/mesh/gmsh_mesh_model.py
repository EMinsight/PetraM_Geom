from __future__ import print_function

import tempfile
import os
import subprocess
import tempfile
import weakref
import numpy as np
import traceback
import time

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GeomModel')

from petram.mesh.mesh_model import Mesh
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

debug = True
class GMesh(Mesh):
    def onItemSelChanged(self, evt):
        '''
        GUI response when model object is selected in
        the dlg_edit_model
        '''
        geom_root = self.geom_root
        if geom_root.is_finalized:
            if geom_root.geom_timestamp != self.geom_timestamp:
                self.onClearMesh(evt)
                self.geom_timestamp = geom_root.geom_timestamp
                evt.Skip()
                return
        viewer = evt.GetEventObject().GetTopLevelParent().GetParent()
        viewer.set_view_mode('mesh', self)
        
    @property
    def geom_timestamp(self):
        return self.parent.geom_timestamp
    
    @geom_timestamp.setter
    def geom_timestamp(self, value):    
        self.parent.geom_timestamp = value
        
        
class GMeshTop(Mesh):
    def attribute_set(self, v):
        v = super(GMeshTop, self).attribute_set(v)
        v['geom_timestamp'] = -1        
        return v
    
class GmshMeshActionBase(GMesh, Vtable_mixin):
    hide_ns_menu = True
    has_2nd_panel = False
    isGmshMesh = True
    
    def attribute_set(self, v):
        v = super(GmshMeshActionBase, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v

    @property
    def geom_root(self):
        return self.root()['Geometry'][self.parent.geom_group]
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

        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
            
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
        
    def onClearMesh(self, evt):
        self.parent.onClearMesh(evt)
        
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
            
    def get_embed(self):
        return [], [], []
            

data = (('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size(def)',
                                default = 1.0, 
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size(def)',
                                default = 1.0, 
                                tip = "CharacteristicLengthMin" )),)
        
                
class GmshMesh(GMeshTop, Vtable_mixin):
    has_2nd_panel = False
    isMeshGroup = True
    vt = Vtable(data)    
#    def __init__(self, *args, **kwargs):
#        super(GmshMesh, self).__init__(*args, **kwargs)
#        NS_mixin.__init__(self, *args, **kwargs)
    @property
    def geom_root(self):
        return self.root()['Geometry'][self.geom_group]
        
    def attribute_set(self, v):
        v['geom_group'] = 'GmshGeom1'
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
        ll.append([None, None, 341, {"label": "Use default size",
                                  "func": 'onSetDefSize',
                                   "noexpand": True}])
        ll.append([None, None, 341, {"label": "Finalize Mesh",
                                  "func": 'onBuildAll',
                                   "noexpand": True}])
        return ll
    
    def get_panel1_value(self):
        return ([self.geom_group,] + list(self.vt.get_panel_value(self)) +
                [self.algorithm, self.algorithm3d, self, self])
    
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
        from .gmsh_mesh_actions import TransfiniteLine, TransfiniteSurface, FreeFace, FreeVolume, FreeEdge, CharacteristicLength, Rotate, Translate, CopyFace, RecombineSurface
        return [FreeVolume, FreeFace, FreeEdge, TransfiniteLine, TransfiniteSurface, CharacteristicLength, Rotate, Translate, CopyFace, RecombineSurface]

    def get_special_menu(self):
        return [('Build All', self.onBuildAll),
                ('Export .geo', self.onExportGeom),
                ('Export .msh', self.onExportMsh),
                ('Clear Mesh', self.onClearMesh)]        
    
    def onSetDefSize(self, evt):
        geom_root = self.geom_root
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
        
        from petram.mesh.gmsh_mesher import write_physical, write_embed
        embed = self.gather_embed()
        geo_text = write_embed(self._txt_rolled[:], embed)        
        geo_text = write_physical(geo_text)
        geo_text.append('Show "*";')
        
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

        
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        geo_text = self._txt_rolled[:] if geo_text is None else geo_text
        
        ret =  generate_mesh(geo_object = None,
                             dim = self._max_mdim,
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text,
                             filename=filename, bin=bin)
        
        viewer.set_figure_data('mesh', self.name(), ret)
        viewer.update_figure('mesh', self.figure_data_name())        
        geom_root = self.geom_root
        viewer._s_v_loop['mesh'] = read_loops(geom_root._txt_unrolled)
        viewer._s_v_loop['geom'] = viewer._s_v_loop['mesh']

        
    def onClearMesh(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        geom_root = self.geom_root

        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
            
        try:
            self.build_mesh(geom_root, nochild =True)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate meshing script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.onUpdateMeshView(evt)
        viewer._view_mode_group = ''
        viewer.set_view_mode('mesh', self)        
        
    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        
        geom_root = self.geom_root
        if not geom_root.is_finalized:
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

        filename = os.path.join(viewer.model.owndir(), self.name())+'_raw'        
        self.onUpdateMeshView(evt, bin='',
                              geo_text = self._txt_rolled[:],
                              filename = filename)
        self.onGenerateMsh(evt)
        evt.Skip()
        
    def onGenerateMsh(self, evt):
        from petram.geom.gmsh_geom_model import read_loops, generate_mesh
        
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        
        filename = os.path.join(viewer.model.owndir(), self.name())
        
        from petram.mesh.gmsh_mesher import write_physical, write_embed
        embed = self.gather_embed()
        geo_text = write_embed(self._txt_rolled[:], embed)        
        geo_text = write_physical(geo_text)

        print("Generating msh with physcal index")
        ret =  generate_mesh(geo_object = None,
                             dim = self._max_mdim,
                             num_quad_lloyd_steps=0,
                             num_lloyd_steps=0,                             
                             geo_text = geo_text,
                             filename=filename, bin='', verbosity='3')


    def gather_embed(self):
        children = [x for x in self.walk()]
        children = children[1:]

        embed = [[], [], []]
        for child in children:
            s, l, p = child.get_embed()
            embed[0].extend(s)
            embed[1].extend(l)
            embed[2].extend(p)
        print("embed", embed)
        return embed
        
    def build_mesh(self, geom_root, stop1=None, stop2=None, filename = None,
                         nochild = False):
        
        from petram.geom.gmsh_geom_model import use_gmsh_api
        
        self.vt.preprocess_params(self)

        if use_gmsh_api:
            from .gmsh_mesh_wrapper import GmshMesher
            mesher = GmshMesher(geom_root,
                            CharacteristicLengthMax = self.clmax,
                            CharacteristicLengthMin = self.clmin,
                            MeshAlgorithm = self.algorithm,
                            MeshAlgorithm3D = self.algorithm3d)                

        else:
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

        if not nochild:
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
        
    def load_gui_figure_data(self, viewer):
        import meshio
        
        filename = os.path.join(viewer.model.owndir(), self.name())+'_raw'
        msh_filename = filename + '.msh'
        if os.path.exists(msh_filename):
            ret = meshio.read(msh_filename)
            return 'mesh', self.name(), ret
        else:
            return 'mesh', self.name(), None
        
    def is_viewmode_grouphead(self):
        return True
    
    def figure_data_name(self):
        try:
            geom_root = self.geom_root
        except:
            return
        if geom_root.is_finalized:
            return self.name(), self.geom_group.strip()
        else:
            print("Geometry not finalized")
            return '', self.geom_group.strip()

    def get_meshfile_path(self):
        '''
        '''
        path = os.path.join(self.root().get_root_path(), self.name() + '.msh')
        if os.path.exists(path):
            dprint1("gmsh file path", path)
            return path
        
        path = os.path.abspath(self.name()+ '.msh')
        if os.path.exists(path):
            dprint1("gmsh file path", path)            
            return path

        assert False, "Mesh file does not exist : " + path

        
    
    
