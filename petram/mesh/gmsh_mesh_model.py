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
        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
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
    
    def onItemSelChanged(self, evt):
        '''
        GUI response when model object is selected in
        the dlg_edit_model
        '''
        geom_root = self.geom_root
        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
        if geom_root.is_finalized:
            if geom_root.geom_timestamp != self.geom_timestamp:
                self.onClearMesh(evt)
                self.geom_timestamp = geom_root.geom_timestamp
                evt.Skip()
                return
        viewer = evt.GetEventObject().GetTopLevelParent().GetParent()
        viewer.set_view_mode('mesh', self)
        
    def get_default_ns(self):
        '''
        this method is overwriten when model wants to
        set its own default namespace. For example, when
        RF module set freq and omega
        '''
        return {"remaining":"remaining", "all":"all"}
    
    
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

        do_clear = True
        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
            
        try:
            filename = os.path.join(viewer.model.owndir(), self.name())+'.msh'            
            kwargs['gui_parent'] = dlg
            kwargs['filename'] = filename
            
            count = self.parent.build_mesh(geom_root, **kwargs)
            do_clear = (count == 0)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate meshing script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.parent.update_meshview(dlg, viewer, clear=do_clear)        
        
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
        # this is default..will be overwitten.
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
        if mode == 'volume':
            viewer.highlight_domain(sel["volume"])
            viewer._dom_bdr_sel = (sel["volume"], [], [], [])
            status_txt = 'Volume :'+ ','.join([str(x) for x in sel["volume"]])
            viewer.set_status_text(status_txt, timeout = 60000)
            viewer._sel_mode = 'volume'
        else:
            figobjs = viewer.highlight_element(sel)
            viewer.set_sel_mode(mode)
            viewer.set_sel_mode() # update buttons        
            if len(figobjs) > 0:
                import ifigure.events
                sel = [weakref.ref(x._artists[0]) for x in figobjs]
                ifigure.events.SendSelectionEvent(figobjs[0], dlg, sel, multi_figobj=figobjs)
            
    def get_embed(self):
        return [], [], []

    def _eval_enitity_id(self, text):
        '''
        "remaining" -> "remaining"
        "all" -> "all"
        something else  -> global vaiable

        failure -> pass thorough
        '''
        if len(text.strip()) == 0: return ''
        try:
             # try to inteprete as integer numbers naively...
             values = list(set([int(x) for x in text.split(',')]))
             values = ','.join([str(x) for x in values])
             return values
        except:
             pass

        # then convert it using namespace                   
        g, l = self.namespace
        ll = {}
        try:                       
            values = eval(text, g, ll)
        except:
            assert False, "can not interpret entity number : " + text                  

        if not isinstance(values, str):
            assert False, "entity id field must be text"

        return values

    def eval_enitity_id(self, *text):
        if len(text) == 1:
            return self._eval_enitity_id(text[0])
        
        return [self._eval_enitity_id(x) for x in text]
    
    
data = (('clmax', VtableElement('clmax', type='float',
                                guilabel = 'Max size(def)',
                                default_txt = '',
                                default = 1e20,
                                tip = "CharacteristicLengthMax" )),
        ('clmin', VtableElement('clmin', type='float',
                                guilabel = 'Min size(def)',
                                default_txt = '',                                
                                default = 0.0, 
                                tip = "CharacteristicLengthMin" )),)
        
                
class GmshMesh(GMeshTop, Vtable_mixin):
    has_2nd_panel = False
    isMeshGroup = True
    vt = Vtable(data)    

    @property
    def geom_root(self):
        return self.root()['Geometry'][self.geom_group]
    @property
    def mesher_data(self):
        if hasattr(self, "_mesher_data"):
            return self._mesher_data
        return None
    
    def attribute_set(self, v):
        v['geom_group'] = 'GmshGeom1'
        v['algorithm'] = 'default'
        v['algorithm3d'] = 'default'
        v['gen_all_phys_entity'] = False
        v['use_profiler'] = False
        v['use_expert_mode'] = False
        super(GmshMesh, self).attribute_set(v)
        self.vt.attribute_set(v)
        return v
    
    def panel1_param(self):    
        ll =   [["Geometry", self.geom_group,  0, {},],]
        ll.extend(self.vt.panel_param(self))

        from petram.mesh.gmsh_mesh_wrapper import Algorithm2D, Algorithm3D

        c1 = list(Algorithm2D)
        c2 = list(Algorithm3D)

        from wx import CB_READONLY
        setting1 = {"style":CB_READONLY, "choices": c1}
        setting2  ={"style":CB_READONLY, "choices": c2}
        
        ll.extend([["2D Algorithm", c1[-1], 4, setting1],
                   ["3D Algorithm", c2[-1], 4, setting2],
                   [None, self.gen_all_phys_entity==1 ,  3,
                    {"text":"Write physical entities for all dimensions."}],
                   [None, self.use_profiler,  3, {"text":"use profiler"}],
                   [None, self.use_expert_mode,  3, {"text":"use GMSH expert mode"}],                   
                   [None, None, 341, {"label": "Use default size",
                                      "func": 'onSetDefSize',
                                      "noexpand": True}], 
                   [None, None, 341, {"label": "Finalize Mesh",
                                      "func": 'onBuildAll',
                                      "noexpand": True}],])


        return ll
    
    def get_panel1_value(self):
        return ([self.geom_group,] + list(self.vt.get_panel_value(self)) +
                [self.algorithm, self.algorithm3d, self.gen_all_phys_entity,
                 self.use_profiler, self.use_expert_mode, self, self, ])
    
    def preprocess_params(self, engine):
        self.vt.preprocess_params(self)
        return
    
    def import_panel1_value(self, v):
        self.geom_group = str(v[0])
        self.vt.import_panel_value(self, v[1:-7])

        self.algorithm = str(v[-7])
        self.algorithm3d = str(v[-6])
        self.gen_all_phys_entity = v[-5]
        self.use_profiler = bool(v[-4])
        self.use_expert_mode = bool(v[-3])        
        
    def panel1_tip(self):
        return ([None] +
                self.vt.panel_tip() +
                ["Alogirth for 2D mesh",
                 "Algoirthm for 3D mesh",
                 "Write lower dimensional physical entity. This may take a long time",
                 "Use cProfiler",
                 "Enable GMSH expert mode to suppress some warning",
                 None, None])
        
    def get_possible_child(self):
        from .gmsh_mesh_actions import TransfiniteLine, TransfiniteSurface, FreeFace, FreeVolume, FreeEdge, CharacteristicLength, CopyFace, CopyFaceRotate, RecombineSurface, ExtrudeMesh, RevolveMesh, MergeText, CompoundCurve, CompoundSurface
        return [FreeVolume, FreeFace, FreeEdge, TransfiniteLine, TransfiniteSurface, CharacteristicLength,  CopyFace, CopyFaceRotate, RecombineSurface, ExtrudeMesh,  RevolveMesh, CompoundCurve, CompoundSurface, MergeText]

    def get_special_menu(self):
        from petram.geom.gmsh_geom_model import use_gmsh_api
        if use_gmsh_api:
            return [('Build All', self.onBuildAll, None),
                    ('Export .msh', self.onExportMsh, None),
                    ('Clear Mesh', self.onClearMesh, None),
                    ('Clear Mesh Sequense...', self.onClearMeshSq, None)]
        else:
             return [('Build All', self.onBuildAll, None),
#                     ('Export .geo', self.onExportGeom, None),
                     ('Export .msh', self.onExportMsh, None),
                     ('Clear Mesh', self.onClearMesh, None),
                     ('Clear Mesh Sequense...', self.onClearMeshSq, None)] 
    
    def onSetDefSize(self, evt):
        geom_root = self.geom_root
        clmax_root, clmin_root = geom_root._clmax_guess
        self.clmax_txt = str(clmax_root)
        self.clmin_txt = str(clmin_root)        
        dlg = evt.GetEventObject().GetTopLevelParent()        
        dlg.OnItemSelChanged()

    '''
    def onExportGeom(self, evt):
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
    '''        
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
        

    def update_meshview(self, dlg, viewer, clear=False):
        import gmsh
        from petram.geom.read_gmsh import read_pts_groups, read_loops

        if clear:
            viewer.del_figure_data('mesh', self.name())
        elif self.mesher_data is None:
            viewer.del_figure_data('mesh', self.name())            
        else:
            print("mesh done face/line",  self._mesh_fface, self._mesh_fline)
            viewer.set_figure_data('mesh', self.name(), self.mesher_data)                

        if 'geom' in viewer._s_v_loop:
            viewer._s_v_loop['mesh'] = viewer._s_v_loop['geom']
        
        viewer.update_figure('mesh', self.figure_data_name())                

        
    def onClearMesh(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        geom_root = self.geom_root

        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)

        do_clear = True
        try:
            count = self.build_mesh(geom_root, nochild =True,
                                    gui_parent = dlg)
            do_clear = (count == 0)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate meshing script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        self.update_meshview(dlg, viewer, clear=do_clear)
        
        viewer._view_mode_group = ''
        viewer.set_view_mode('mesh', self)
        
    def onClearMeshSq(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        
        import ifigure.widgets.dialog as dialog
        ret = dialog.message(parent=dlg,
                             message='Are you sure to delete all Mesh sequence',
                             title='Mesh Sequence Delete',
                             style=2)
        if ret == 'ok':
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate meshing script',
                                 title='Error',
                                 traceback=traceback.format_exc())
            
            for x in self.keys():
                del self[x]
            dlg.tree.RefreshItems()

    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.import_selected_panel_value()
        
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        
        geom_root = self.geom_root
        if not geom_root.is_finalized:
            geom_root.onBuildAll(evt)
            
        do_clear = True            
        try:
            filename = os.path.join(viewer.model.owndir(), self.name())+'.msh'
            count = self.build_mesh(geom_root, finalize=True, filename=filename,
                                    gui_parent = dlg)
            do_clear = count == 0
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to generate mesh script',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()


        self.update_meshview(dlg, viewer, clear=do_clear)
        
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

    def build_mesh(self, geom_root, stop1=None, stop2=None, filename = '', 
                         nochild = False, finalize=False, gui_parent = None):
        import gmsh
        from petram.geom.read_gmsh import read_pts_groups, read_loops
        
        self.vt.preprocess_params(self)
        clmax, clmin = self.vt.make_value_or_expression(self)        
        dprint1("calling build mesh with", clmax, clmin)
        geom_root = self.geom_root
        
        if not geom_root.is_finalized:
            geom_root.build_geom4(no_mesh=True, finalize=True)
            
        #else:
        #    geom = geom_root._gmsh4_data[-1]

        from petram.mesh.gmsh_mesh_wrapper import GMSHMeshWrapper as GmshMesher
        mesher = GmshMesher(format=2.2,
                            CharacteristicLengthMax = clmax,
                            CharacteristicLengthMin = clmin,
                            EdgeResolution = 3,                             
                            MeshAlgorithm = self.algorithm,
                            MeshAlgorithm3D = self.algorithm3d,
                            use_profiler = self.use_profiler,
                            use_expert_mode = self.use_expert_mode,
                            gen_all_phys_entity = self.gen_all_phys_entity)
        
        #mesher.load_brep(geom_root._geom_brep)
        
        children = [x for x in self.walk()]
        children = children[1:]

        if not nochild:
            for child in children:
                if not child.enabled: continue
                if child is stop1: break            # for build before                
                child.vt.preprocess_params(child)
                #child.check_master_slave(mesher)                
                child.add_meshcommand(mesher)
                if child is stop2: break            # for build after

        if mesher.count_sequence() > 0:
            self._mesher_data = None                # set None since mesher may die...
            
            L = mesher.count_sequence()*4 + 3
            
            if gui_parent is not None:
                import wx                
                gui_parent = wx.GetApp().TopWindow
                pgb = wx.ProgressDialog("Generating mesh...",
                                "", L, parent = gui_parent,
                                style = wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
                def close_dlg(evt, dlg=pgb):
                    pgb.Destroy()
                pgb.Bind(wx.EVT_CLOSE, close_dlg)
            else:
                pgb = None
            
            max_mdim, done, data, msh_output = mesher.run_generater(geom_root._geom_brep,
                                                                    filename, 
                                                                    finalize=finalize,
                                                                    progressbar = pgb)            
            if pgb is not None: pgb.Destroy()
            
            self._mesher_data = data         
            self._max_mdim = max_mdim
            if finalize:
                self._msh_output = msh_output
            else:
                self._msh_output = ''
        else:
            self._max_mdim = 0
            done = [], [], [], []
            
        self._mesh_fface = done[2] # finished surfaces
        self._mesh_fline = done[1] # finished lines

        
        return (mesher.count_sequence() > 0)
    
    def generate_mesh_file(self):
        cwd = os.getcwd()
        dprint1("Generating Mesh in " + cwd)
        geom_root = self.geom_root
        filename = os.path.join(cwd, self.name())+'.msh'
        count = self.build_mesh(geom_root, finalize=True, filename=filename,
                                gui_parent = None)
        if count == 0:
            assert False, "Failed to generate mesh"
        else:
            dprint1("Generating Mesh ... Done")        

              
    def load_gui_figure_data(self, viewer):
        #import meshio
        return 'mesh', self.name(), None
    
        #filename = os.path.join(viewer.model.owndir(), self.name())+'_raw'
        #msh_filename = filename + '.msh'
        #if os.path.exists(msh_filename):
        #    ret = meshio.read(msh_filename)
        #    return 'mesh', self.name(), ret
        #else:

        
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
        if hasattr(self, '_msh_output') and self._msh_output != '':
            path = self._msh_output
            if os.path.exists(path):
                dprint1("gmsh file path", path)
                return path
                
        path = os.path.join(self.root().get_root_path(), self.name() + '.msh')
        if os.path.exists(path):
            dprint1("gmsh file path", path)
            return path
        
        path = os.path.abspath(self.name()+ '.msh')
        if os.path.exists(path):
            dprint1("gmsh file path", path)            
            return path

        assert False, "Mesh file does not exist : " + path

        
    
    
