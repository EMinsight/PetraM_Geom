'''

    Model Tree to stroe MFEM model parameters

'''
from __future__ import print_function

import tempfile
import os
import subprocess
import traceback
import sys
import re
import time
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import warnings

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshGeomModel')

from petram.model import Model
import time

from threading import Thread
try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty  # python 3.x

import petram.geom.gmsh_config
from petram.geom.geom_model import GeomBase, GeomTopBase
from petram.namespace_mixin import NS_mixin
from petram.phys.vtable import VtableElement, Vtable, Vtable_mixin

debug = True
geom_key_dict = {'SurfaceBase': 'sb',
                 'PlaneSurface' : 'sp',
                 'Point': 'pt',
                 'Line': 'ln',
                 'Spline': 'sp',
                 'SurfaceID': 's',
                 'VolumeID': 'v',
                 'LineID': 'l',
                 'VertexID': 'p'}


def get_gmsh_exe():
    macos_gmsh_location = '/Applications/Gmsh.app/Contents/MacOS/gmsh'
    if os.path.isfile(macos_gmsh_location):
        gmsh_executable = macos_gmsh_location
    else:
        gmsh_executable = 'gmsh'
    return gmsh_executable


def get_gmsh_major_version():
    gmsh_exe = get_gmsh_exe()
    try:
        out = subprocess.check_output(
                [gmsh_exe, '--version'],
                 stderr=subprocess.STDOUT
                 ).strip().decode('utf8')
    except:
        return -1
    ex = out.split('.')
    return int(ex[0])

use_gmsh_api = True
gmsh_Major=get_gmsh_major_version()
if gmsh_Major <= 3: use_gmsh_api = False

def enqueue_output(p, queue):
    while True:
        line = p.stdout.readline()
        
        queue.put(line.strip())
        if p.poll() is not None: 
            queue.put("End of Thread")
            return
        print(line.strip())
    queue.put("End of Thread")    
    
def collect_std_out(p,  verbose=True):
    q = Queue()
    t = Thread(target=enqueue_output, args=(p, q))
    t.daemon = True # thread dies with the program
    t.start()
    
    lines = []
    alive = True
    while True:

        time.sleep(0.01)
        
        try:  line = q.get_nowait() # or q.get(timeout=.1)
        except Empty:
            if p.poll() is not None:
                print('proces exited')
                break
            else:
                continue
        ec = p.poll()
        if ec is not None and ec < 0:
            print("RETURNIng due to this?")
            break  # on unix, this means process killed by a signal
        lines.append(line)
    return lines, p.poll()

def get_geom_key(obj):
    if obj.__class__.__name__ in geom_key_dict:
        return geom_key_dict[obj.__class__.__name__]

    assert False, " name not found for " + obj.__class__.__name__ + " in " + str(geom_key_dict)

    ### it should not come here...
    name = obj.__class__.__name__    
    key = ''.join([i.lower() for i in name if not i.isupper()])
    for k in geom_key_dict.keys():
        if geom_key_dict[k] == key:
            assert False, key + " is used for " + k.__name__
            
    geom_key_dict[obj.__class__] = key
    if debug: print(geom_key_dict)
    return key

def get_twoletter_keys(t):
    if t == 'p': return 'pt'
    elif t == 'l': return 'ln'    
    elif t == 'f': return 'fs'
    elif t == 'v': return 'vl'
    else:
        return t

class GeomObjs(dict):
    def duplicate(self):
        if not hasattr(self, "_past_keys"):
            self._past_keys = []            
        obj = GeomObjs(self)
        obj._past_keys = self._past_keys
        return obj
        
    def addobj(self, obj, name):
        key = ''.join([i for i in name if not i.isdigit()])
        key = get_twoletter_keys(key)        
        if not hasattr(self, "_past_keys"):
            self._past_keys = []
        keys = self._past_keys
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
        self._past_keys.append(newkey)
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

    #def build_geom(self, geom, objs):
    #    self._newobjs = []
    #    warnings.warn("Not implemented: " + self.__class__.__name__, Warning)        

    def gsize_hint(self, geom, objs):
        '''
        return quick estimate of geometry size min and max
        '''
        warnings.warn("Not implemented", Warning)
        return -1, -1
    
    def get_special_menu(self):
        return [('Build this step', self.onBuildAfter, None)]

    def _onBuildThis(self, evt, **kwargs):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()
        kwargs['gui_parent'] = dlg
        try:
            p  = self.parent
            if isinstance(p, GmshGeom):
                rootg = p
            else: # work plane
                rootg = p.parent
            rootg._geom_finalized = False
            
            od = os.getcwd()
            os.chdir(viewer.model.owndir())
            rootg.build_geom(**kwargs)
            os.chdir(od)            

        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to build geometry',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()
        rootg.onUpdateGeoView(evt)

        
    def onBuildBefore(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.import_selected_panel_value()
        
        mm = dlg.get_selected_mm()
        if mm is None: return

        self._onBuildThis(evt, stop1 = mm)
        evt.Skip()
        
    def onBuildAfter(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.import_selected_panel_value()
        mm = dlg.get_selected_mm()
        if mm is None: return
        
        self._onBuildThis(evt, stop2 = mm)

        dlg = evt.GetEventObject().GetTopLevelParent()
        dlg.select_next_enabled()
        evt.Skip()

    def add_geom_sequence(self, geom):
        gui_name = self.fullname()
        self.vt.preprocess_params(self)                
        gui_param = self.vt.make_value_or_expression(self)
        geom_name = self.__class__.__name__
        geom.add_sequence(gui_name, gui_param, geom_name)
'''        
class BrepFile(GeomTopBase):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(BrepFile, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    @property
    def is_finalized(self):
        if not hasattr(self, "_geom_finalized"):
            self._geom_finalized = False
        return self._geom_finalized
    
    @property
    def _geom_brep(self):
        return self.brep_file_path

    @property
    def geom_finalized(self):
        if not hasattr(self, "_geom_finalized"):
            self._geom_finalized = False
        return self._geom_finalized
    
    @geom_finalized.setter
    def geom_finalized(self, value):
        self._geom_finalized = value
    
    def attribute_set(self, v):
        v = super(BrepFile, self).attribute_set(v)
        v['brep_file_path'] = ''
        v['geom_timestamp'] = 0
        v['geom_prev_algorithm'] = 2
        v['geom_prev_res'] = 30
        return v

    def panel1_param(self):
        import wx
        
        wc = "ANY|*|Brep|*.brep"
        return [["File(.brep)", None, 45, {'wildcard':wc}],        
                ["PreviewAlgorith", "Automatic", 4, {"style":wx.CB_READONLY,
                                                     "choices": ["Auto", "MeshAdpat",
                                                                 "Delaunay", "Frontal"]}],
                ["PreviewResolution", 30,  400, None],
                [None, None, 341, {"label": "Finalize Geom",
                                   "func": 'onBuildAll',
                                   "noexpand": True}],]

    def get_panel1_value(self):
        aname = {2: "Auto", 1: "MeshAdpat", 5: "Delaunay", 6:"Frontal"}
        txt = aname[self.geom_prev_algorithm]
        return [self.brep_file_path, txt, self.geom_prev_res, self]
       
    def import_panel1_value(self, v):
        aname = {2: "Auto", 1: "MeshAdpat", 5: "Delaunay", 6:"Frontal"}
        for k in aname:
            if v[1] == aname[k]:
                self.geom_prev_algorithm = k

        self.geom_prev_res = int(v[2])
        self.brep_file_path = str(v[0])
        
    def get_special_menu(self):
        return [('Load File', self.onBuildAll, None),]
    
    def onBuildAll(self, evt):
        import gmsh
        
        if not hasattr(self, "_gmsh4_data"):
            self._gmsh4_data = None
        #if self._gmsh4_data is not  None:
        #    self._gmsh4_data[-1].finalize()
            
        objs = GeomObjs()
        self._objs = objs

        from petram.geom.gmsh_geom_wrapper import Geometry
        geom = Geometry()
        
        import gmsh
        geom.clear()        
        gmsh.model.occ.importShapes(self.brep_file_path, highestDimOnly=False)
        gmsh.model.occ.synchronize()

            
        # here we ask for 2D mesh for plotting.
        # size control is done based on geometry size.
        size = []
        for dim, tag in gmsh.model.getEntities():
            x1, y1, z1, x2, y2, z2 = gmsh.model.getBoundingBox(dim, tag)
            s = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
            size.append((dim, tag, s))
        maxsize = max([x[-1] for x in size])
        lcar = defaultdict(lambda: maxsize)
        for dim, tag in gmsh.model.getEntities(1):
            x1, y1, z1, x2, y2, z2 = gmsh.model.getBoundingBox(dim, tag)
            s = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
            bdimtags = gmsh.model.getBoundary(((dim, tag,),), oriented=False)
            for bdim, btag in bdimtags:
                lcar[btag] = min((lcar[btag], s))
        lcar = dict(lcar)
        print(lcar)
        #dim2_size = min([s[2] for s in ss if s[0]==2]+[3e20])
        #dim1_size = min([s[2] for s in ss if s[0]==1]+[3e20])

        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", maxsize/self.geom_prev_res)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min(lcar))
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1)
        geom.model.mesh.generate(1)
        
        gmsh.option.setNumber("Mesh.Algorithm", self.geom_prev_algorithm)
        #gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1e22)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",maxsize/self.geom_prev_res)
        gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)        
        geom.model.mesh.generate(2)        


        self.geom_finalized = True
        self.geom_timestamp = time.ctime()        
        
        from petram.geom.read_gmsh import read_pts_groups, read_loops
        ptx, cells, cell_data = read_pts_groups(geom)
        l, s, v = read_loops(gmsh)
        self._gmsh4_data = (ptx, cells, cell_data, l, s, v, geom)
        
        ret = ptx, cells, {}, cell_data, {}

        # set clmax guess from geometry size
        clmax = max(lcar)/3.
        clmin = min(lcar)/3.
        self._clmax_guess = (clmax, clmin)
        self._geom_coords = ret

        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        
        viewer.set_figure_data('geom', self.name(), ret)
        viewer.update_figure('geom', self.name())
        
        viewer._s_v_loop['geom'] = s, v
        viewer._s_v_loop['mesh'] = s, v

        evt.Skip()
'''        

class GmshGeom(GeomTopBase):
    has_2nd_panel = False
    def __init__(self, *args, **kwargs):
        super(GmshGeom, self).__init__(*args, **kwargs)
        NS_mixin.__init__(self, *args, **kwargs)
        
    @property
    def is_finalized(self):
        if not hasattr(self, "_geom_finalized"):
            self._geom_finalized = False
        return self._geom_finalized
    
    @property
    def geom_finalized(self):
        if not hasattr(self, "_geom_finalized"):
            self._geom_finalized = False
        return self._geom_finalized
    
    @geom_finalized.setter
    def geom_finalized(self, value):
        self._geom_finalized = value
    
    @property    
    def build_stop(self):
        if not hasattr(self, "_build_stop"):
            self._build_stop = (None, None)
        return self._build_stop
    
    def attribute_set(self, v):
        v = super(GmshGeom, self).attribute_set(v)
        v['geom_timestamp'] = 0
        v['geom_prev_algorithm'] = 2
        v['geom_prev_res'] = 3
        v['occ_parallel'] = False
        v['maxthreads'] = 1
        v['skip_final_frag'] = False
        v['use_1d_preview'] = False
        v['use_curvature'] = False
        v['long_edge_thr'] = 0.3        
        v['small_edge_thr'] = 0.001
        v['small_edge_seg'] = 1
        v['max_seg'] = 30
        return v
        
    def get_possible_child(self):
        from .gmsh_primitives import Point, Line, Spline, Circle, Rect, Polygon, Box, Ball, Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep, LineLoop, CreateLine, CreateSurface, CreateVolume, SurfaceLoop, Union, Intersection, Difference, Fragments, Copy, Remove, Move, Rotate, Flip, Scale, WorkPlane, WorkPlaneByPoints, healCAD, CADImport, BrepImport,  Fillet, Chamfer, Array, ArrayRot
        return [Point,  Line, Circle, Rect, Polygon, Spline, Box, Ball, Cone, Wedge, Cylinder, Torus, CreateLine, CreateSurface, CreateVolume, LineLoop, SurfaceLoop, Extrude, Revolve, Sweep, Union, Intersection, Difference, Fragments, Copy, Remove, Move, Rotate, Flip, Scale, WorkPlane, WorkPlaneByPoints, healCAD, CADImport, BrepImport, Fillet, Chamfer, Array, ArrayRot]
    
    def get_possible_child_menu(self):
        from .gmsh_primitives import Point, Line, Spline, Circle, Rect, Polygon, Box, Ball, Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep, LineLoop, CreateLine, CreateSurface, CreateVolume, SurfaceLoop, Union, Intersection, Difference, Fragments, Copy, Remove, Move, Rotate, Flip, Scale, WorkPlane, WorkPlaneByPoints, healCAD, CADImport, BrepImport, Fillet, Chamfer, Array, ArrayRot
        return [("", Point),("", Line), ("", Circle), ("", Rect), ("", Polygon),
                ("", Spline),("", Fillet), ("", Chamfer), 
                ("3D shape...", Box),
                ("", Ball), ("", Cone), ("", Wedge), ("", Cylinder),
                ("!", Torus),
                ("", CreateLine), ("", CreateSurface), ("", CreateVolume),
                #("", LineLoop), ("", SurfaceLoop),
                ("Protrude...", Extrude), ("", Revolve), ("!", Sweep),
                ("", Copy), ("", Remove),
                ("Translate...", Move,), ("", Rotate),("", Flip),("", Scale),
                ("", Array), ("!", ArrayRot),
                ("Boolean...", Union),("",Intersection),("",Difference),("!",Fragments),
                ("WorkPlane...", WorkPlane), ("!", WorkPlaneByPoints),
                ("Import...", BrepImport),("", CADImport),("!", healCAD),
                ]
                
    def get_special_menu(self):
        if use_gmsh_api:
            return [('Build All', self.onBuildAll, None),
                    ('Export .brep', self.onExportBrep, None)]
        else:
            return [('Build All', self.onBuildAll, None),
                    ('Export .geo', self.onExportGeom, None)]
    
    def panel1_param(self):
        import wx
        return [["", "Geometry model using GMSH", 2, None],
                ["PreviewAlgorith", "Automatic", 4, {"style":wx.CB_READONLY,
                                                     "choices": ["Auto", "MeshAdpat",
                                                                 "Delaunay", "Frontal"]}],
                ["Preview Resolution", 30,  400, None],
                ["Long  Edge Thr.", self.long_edge_thr, 300, None],                
                ["Small Edge Thr.", self.small_edge_thr, 300, None],
                ["Small Edge #Seg.", self.small_edge_seg, 400, None],
                ["Max #seg in Preview", self.max_seg, 400, None],
                ["Preview #threads", self.maxthreads, 400, None],
                [None, self.occ_parallel, 3, {"text":"OCC parallel boolean"}],
                [None, self.skip_final_frag, 3, {"text":"Skip fragmentationn"}],
                [None, self.use_1d_preview, 3, {"text":"Use line preview"}],
                [None, self.use_curvature, 3, {"text":"Consider curvature in preview generation"}],                                                
                [None, None, 341, {"label": "Finalize Geom",
                                   "func": 'onBuildAll',
                                   "noexpand": True}],]

    def get_panel1_value(self):
        aname = {2: "Auto", 1: "MeshAdpat", 5: "Delaunay", 6:"Frontal"}
        txt = aname[self.geom_prev_algorithm]
        return [None, txt, self.geom_prev_res, self.long_edge_thr,
                self.small_edge_thr, self.small_edge_seg, self.max_seg,
                self.maxthreads, self.occ_parallel,
                self.skip_final_frag, self.use_1d_preview, self.use_curvature, self]
       
    def import_panel1_value(self, v):
        aname = {2: "Auto", 1: "MeshAdpat", 5: "Delaunay", 6:"Frontal"}
        for k in aname:
            if v[1] == aname[k]:
                self.geom_prev_algorithm = k

        self.geom_prev_res = int(v[2])
        self.long_edge_thr = float(v[3])
        self.small_edge_thr = float(v[4])
        self.small_edge_seg = int(v[5])
        self.max_seg = int(v[6])                        
        self.maxthreads  =  int(v[7])
        self.occ_parallel  = v[8]
        self.skip_final_frag = v[9]
        self.use_1d_preview = v[10]
        self.use_curvature = v[11]

    def onBuildAll(self, evt):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        engine = viewer.engine
        engine.build_ns()

        try:
            od = os.getcwd()
            os.chdir(viewer.model.owndir())
            self.build_geom(finalize = True, gui_parent=dlg)
            os.chdir(od)
        except:
            import ifigure.widgets.dialog as dialog               
            dialog.showtraceback(parent = dlg,
                                 txt='Failed to build geometry',
                                 title='Error',
                                 traceback=traceback.format_exc())
        dlg.OnRefreshTree()

        filename = os.path.join(viewer.model.owndir(), self.name())
        self.onUpdateGeoView(evt, filename = filename)

        '''
        if not use_gmsh_api:
            fid = open(filename + '.geo_unrolled', 'w')
            fid.write('\n'.join(self._txt_unrolled))
            fid.close()
        '''    
        self.geom_finalized = True
        self.geom_timestamp = time.ctime()
        evt.Skip()
        
    def onUpdateGeoView4(self, evt, filename = None):
        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        ptx, cells, cell_data, l, s, v, geom = self._gmsh4_data
        ret = ptx, cells, {}, cell_data, {}

        # set clmax guess from geometry size
        #xmin, ymin, zmin, xmax, ymax, zmax = geom.getBoundingBox()
        #l = ((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)**0.5
        #clmax = l/3.
        #clmin = l/300.
        #self._clmax_guess = (clmax, clmin)
        
        self._geom_coords = ret
        viewer.set_figure_data('geom', self.name(), ret)
        viewer.update_figure('geom', self.name())
        
        viewer._s_v_loop['geom'] = s, v
        viewer._s_v_loop['mesh'] = s, v
        
    def onUpdateGeoView(self, evt, filename = None):       
        if globals()['gmsh_Major']==4 and use_gmsh_api:
            return self.onUpdateGeoView4(evt, filename = filename)
        else:
            assert False, "GMSH 3 is not supported"

    def walk_over_geom_chidlren(self, geom, stop1=None, stop2=None):
        self._build_stop = (None, None)
        
        children = [x for x in self.walk()]
        children = children[1:]
        for child in children:
            if hasattr(child, "_newobjs"): del child._newobjs
            
        children = self.get_children()
        for child in children:
            if not child.enabled: continue
            
            if len(child.get_children())==0:
                child.vt.preprocess_params(child)
                if child is stop1: break            # for build before
                child.add_geom_sequence(geom)
                if child is stop2: break            # for build after
                
            else:  # workplane
                children2 = child.get_children()
                child.vt.preprocess_params(child)
                if child is stop1: break            # for build before                

                do_break = False
                
                geom.add_sequence('WP_Start', 'WP_Start', 'WP_Start')
                for child2 in children2:
                    if not child2.enabled: continue                    
                    child2.vt.preprocess_params(child2)
                    if child2 is stop1:
                        do_break = True
                        break            # for build before
                    child2.add_geom_sequence(geom)                    
                    if child2 is stop2:
                        do_break = True                        
                        break            # for build after

                # translate 2D objects in 3D space
                #for x in org_keys: del objs2[x]
                child.add_geom_sequence(geom)
                geom.add_sequence('WP_End', 'WP_End', 'WP_End')
                
                if do_break: break
                if child is stop2: break            # for build after
        if stop1 is not None:
            self._build_stop = (stop1, None)
            return stop1.name()
        if stop2 is not None:
            self._build_stop = (None, stop2)
            return stop2.name()
        return self.name()

    def update_GUI_after_geom(self, data, objs):
        children = [x for x in self.walk()]
        children = children[1:]
        
        for child in children:
            if hasattr(child, "_newobjs"): del child._newobjs

        for child in children:
            if child.fullname() in data:
                dd = data[child.fullname()]
                child._objkeys = dd[0]
                child._newobjs = dd[1]  
                
        self._objs = objs        

    def check_create_new_child(self,gs):
        if not hasattr(self, '_prev_sequence'):
            return True

        if len(gs.geom_sequence) < len(self._prev_sequence):
            return True

        import six        
        if six.PY2:
            import cPickle as pickle
        else:
            import pickle
        
        for k, s in enumerate(self._prev_sequence):
            s_txt1 = pickle.dumps(s)
            s_txt2 = pickle.dumps(gs.geom_sequence[k])
            if s_txt1 != s_txt2:
                return True
            else:
                dprint1("check passed", s[0])
        return False
        
    def build_geom4(self, stop1=None, stop2=None, filename = None,
                    finalize = False, no_mesh=False, gui_parent=None):
        '''
        filename : export geometry to a real file (for debug)
        '''
        import gmsh
        
        if not hasattr(self, "_gmsh4_data"):
            self._gmsh4_data = None
        #if self._gmsh4_data is not  None:
        #    self._gmsh4_data[-1].finalize()

        from petram.geom.gmsh_geom_wrapper import GeometrySequence,GMSHGeometryGenerator
        '''
        geom = Geometry(PreviewResolution = self.geom_prev_res,
                        PreviewAlgorithm = self.geom_prev_algorithm,
                        OCCParallel = int(self.occ_parallel),
                        Maxthreads = self.maxthreads,
                        SkipFrag = self.skip_final_frag,
                        Use1DPreview = self.use_1d_preview)

        geom.set_factory('OpenCASCADE')
        '''

        gs = GeometrySequence()
        stopname = self.walk_over_geom_chidlren(gs, stop1=stop1, stop2=stop2)

        import wx
        if gui_parent is None:
            gui_parent = wx.GetApp().TopWindow

        L = len(gs.geom_sequence) + 3
        pgb = wx.ProgressDialog("Generating geometry...",
                                "", L, parent = gui_parent,
                                style = wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT)
        def close_dlg(evt, dlg=pgb):
            pgb.Destroy()
        pgb.Bind(wx.EVT_CLOSE, close_dlg)

        if (hasattr(self, "_p") and self._p[0].is_alive()):
            new_process = self.check_create_new_child(gs)
        else:
            new_process = True

        if new_process:
            if (hasattr(self, "_p") and self._p[0].is_alive()):
                self._p[0].terminate()

            task_q = mp.Queue() # data to child
            q =  mp.Queue() # data from child
            p = GMSHGeometryGenerator(q, task_q)
            p.start()
            print("process ID", p.pid)
            self._p = (p, task_q, q)
            self._prev_sequence = gs.geom_sequence
            start_idx = 0
        else:
            ll = len(self._prev_sequence)
            self._prev_sequence = gs.geom_sequence
            start_idx = ll
            
        success, dataset  = gs.run_generator(self, no_mesh = no_mesh, finalize=finalize,
                                             filename = stopname, progressbar = pgb,
                                             create_process = new_process,
                                             process_param = self._p,
                                             start_idx = start_idx)

        pgb.Destroy()

        if not success:
            print(dataset) # this is an error message
            self._p[0].terminate()
            self._prev_sequence = []
            return
        
        gui_data, objs, brep_file, data, vcl = dataset

        self._geom_brep = brep_file
        self.update_GUI_after_geom(gui_data, objs)

        if data is None:  # if no_mesh = True
            return   
        # for the readablity I expend data here, do we need geom?        
        ptx, cells, cell_data, l, s, v = data
        self._gmsh4_data = (ptx, cells, cell_data, l, s, v, gs)

        self._clmax_guess = vcl

        return


    def build_geom(self, stop1=None, stop2=None, filename = None,
                   finalize = False, gui_parent=None):

        if globals()['gmsh_Major']==4 and use_gmsh_api:
            self.build_geom4(stop1=stop1, stop2=stop2,
                             filename=filename,
                             finalize=finalize,
                             gui_parent = gui_parent)
        else:
            assert False, "GMSH 3 is not supported"
            #self.build_geom3(stop1=stop1, stop2=stop2,
            #                 filename=filename,
            #                 finalize=finalize)

    def onExportGeom(self, evt):
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
            
    def onExportBrep(self, evt):
        if not hasattr(self, "_geom_brep"):
            evt.Skip()
            return
        from ifigure.widgets.dialog import write
        parent = evt.GetEventObject()
        path = write(parent,
                     message = 'Enter .brep file name',
                     wildcard = '*.brep')
        if path != '':
            from shutil import copyfile
            copyfile(self._geom_brep, path)

    '''
    def load_gui_figure_data(self, viewer):
        import meshio
        filename = os.path.join(viewer.model.owndir(), self.name())
        msh_filename = filename + '.msh'
        if not os.path.exists(msh_filename):
            return 'geom', self.name(), None
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
    '''
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

