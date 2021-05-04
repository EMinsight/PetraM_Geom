import numpy as np


from petram.geom.gmsh_geom_model import GmshGeom

class OCCGeom(GmshGeom):
    has_2nd_panel = False

    @classmethod
    def fancy_menu_name(cls):
        return 'OCC Geometry'

    @classmethod
    def fancy_tree_name(cls):
        return 'OCCSequence'

    def attribute_set(self, v):
        v = super(OCCGeom, self).attribute_set(v)
        v['long_edge_thr'] = 0.5
        v['small_edge_thr'] = 0.01
        v['use_occ_preview'] = True
        return v

    def build_geom4(self, stop1=None, stop2=None, filename=None,
                    finalize=False, no_mesh=False, gui_parent=None,
                    cwd=None):

        self.use_occ_preview = True
        self.do_build_geom4(stop1=stop1, stop2=stop2, filename=filename,
                            finalize=finalize, no_mesh=no_mesh, gui_parent=gui_parent,
                            cwd=cwd)

    def inspect_geom(self, inspect_type, params):
        '''
        command = (type_of_inspect, params)
        '''
        if not hasattr(self, '_gso'):
            assert False, "Geometry Sequence Operator does not exist"
        command = (inspect_type, params)
        return self._gso.inspect_geom(command)

    def onExportSelectedBrep(self, evt):
        if not hasattr(self, '_gso'):
            return None

        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()

        selection = viewer.dom_bdr_sel
        kind = viewer.get_sel_mode()
        if kind == 'volume':
            selection = selection[0], [], [], []
        elif kind == 'face':
            selection = [], selection[1], [], []
        elif kind == 'edge':
            selection = [], [], selection[2], []
        elif kind == 'point':
            selection = [], [], [], selection[3]
        else:
            return None

        parent = evt.GetEventObject()
        
        from ifigure.widgets.dialog import write        
        path = write(parent,
                     message='Enter .brep file name',
                     wildcard='*.brep')

        if path != '':
            if not path.endswith('.brep'):
                path = path + '.brep'
            return self._gso.export_shapes(selection, path)

        evt.Skip()
        return None

    def panel1_param(self):
        import wx
        return [["", "Geometry model using OpenCascade", 2, None],
                #["PreviewAlgorith", "Automatic", 4, {"style": wx.CB_READONLY,
                #                                     "choices": ["Auto", "MeshAdpat",
                #                                                 "Delaunay", "Frontal"]}],
                ["Preview Resolution (linear)", self.small_edge_thr, 300, None],
                ["Preview Resolution (angle)", self.long_edge_thr, 300, None],
                [None, self.maxthreads > 1, 3, {"text": "Parallel preview"}],                
                [None, self.occ_parallel, 3, {"text": "Parallel boolean"}],
                [None, self.skip_final_frag, 3, {
                    "text": "Skip fragmentationn"}],
                [None, None, 341, {"label": "Finalize Geom",
                                   "func": 'onBuildAll',
                                   "noexpand": True}], ]

    def get_panel1_value(self):
        return [None, self.small_edge_thr, self.long_edge_thr,
                self.maxthreads, self.occ_parallel, self.skip_final_frag, self]

    def import_panel1_value(self, v):
        self.small_edge_thr = float(v[1])
        self.long_edge_thr = float(v[2])
        self.maxthreads = 2 if v[3] else 1
        self.occ_parallel = v[4]
        self.skip_final_frag = v[5]
        self.use_occ_preview = True

    def get_possible_child(self):
        from petram.geom.geom_primitives import (PointOCC, LineOCC, CircleOCC, Polygon2,
                                                 Point, PointCenter, PointByUV, PointOnEdge,
                                                 PointCircleCenter, Line, Spline,
                                                 Circle, CircleByAxisPoint, CircleBy3Points,
                                                 CircleByAxisCenterRadius,
                                                 Rect, Polygon, OCCPolygon, Box, Ball,
                                                 Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep,
                                                 LineLoop, CreateLine, CreateSurface, CreateVolume,
                                                 SurfaceLoop, Union, Union2, MergeFace,
                                                 Intersection, Difference, Fragments,
                                                 SplitByPlane, Copy, Remove, Remove, Remove2, RemoveFaces,
                                                 Move, Rotate, Flip, Scale, WorkPlane,
                                                 WorkPlaneByPoints, WPParallelToPlane,
                                                 WPNormalToPlane,
                                                 healCAD, CADImport, BrepImport,
                                                 Fillet, Chamfer,
                                                 Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                                                 ArrayPath,
                                                 ThruSection, CreateShell,
                                                 RotateCenterPoints, MoveByPoints, ExtendedLine)
        return [PointOCC, LineOCC, CircleOCC, Polygon2,
                Point, PointCenter, PointOnEdge, PointByUV, PointCircleCenter,
                Line, Circle, CircleByAxisPoint, CircleBy3Points,
                CircleByAxisCenterRadius,
                Rect, Polygon,  OCCPolygon, Spline, Box,
                Ball, Cone, Wedge, Cylinder,
                Torus, CreateLine, CreateSurface, CreateVolume, LineLoop, SurfaceLoop,
                Extrude, Revolve, Sweep, Union, Union2, MergeFace,
                Intersection, Difference, Fragments, SplitByPlane, Copy, Remove,
                Remove2, RemoveFaces, Move, Rotate,
                Flip, Scale, WorkPlane, WorkPlaneByPoints, WPParallelToPlane,
                WPNormalToPlane,
                healCAD, CADImport, BrepImport,
                Fillet, Chamfer, Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                ArrayPath,     
                ThruSection, CreateShell, RotateCenterPoints, MoveByPoints, ExtendedLine]

    def get_possible_child_menu(self):
        from petram.geom.geom_primitives import (PointOCC, LineOCC, CircleOCC, Polygon2,
                                                 Point, PointCenter, PointCircleCenter,
                                                 PointOnEdge, PointByUV, Line, Spline,
                                                 Circle, CircleByAxisPoint, CircleBy3Points,
                                                 CircleByAxisCenterRadius,
                                                 Rect, Polygon, OCCPolygon, Box, Ball,
                                                 Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep,
                                                 LineLoop, CreateLine, CreateSurface, CreateVolume,
                                                 SurfaceLoop, Union, Union2, MergeFace,
                                                 Intersection, Difference, Fragments,
                                                 SplitByPlane, Copy, Remove, Remove2, RemoveFaces,
                                                 Move, Rotate, Flip, Scale,
                                                 WorkPlane, WorkPlaneByPoints, WPParallelToPlane,
                                                 WPNormalToPlane,
                                                 healCAD, CADImport, BrepImport,
                                                 Fillet, Chamfer,
                                                 Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                                                 ArrayPath,
                                                 ThruSection, CreateShell,
                                                 RotateCenterPoints, MoveByPoints, ExtendedLine)

        return [("Geometry Element...", None),
                ("Points...", PointOCC), ("", PointCenter), ("", PointOnEdge),
                ("", PointCircleCenter), ("!", PointByUV),
                ("Lines...", LineOCC), ("!", ExtendedLine),
                ("Polygon...", Polygon2), ("!", OCCPolygon), 
                ("Circle...", CircleOCC), ("", CircleByAxisPoint),
                ("", CircleByAxisCenterRadius), ("!", CircleBy3Points),
                ("", Rect),
                ("", Spline),
                ("!", None),
                ("3D shape...", Box),
                ("", Ball), ("", Cone), ("", Wedge), ("", Cylinder),
                ("!", Torus),
                ("Create...", CreateLine), ("", CreateSurface), ("", CreateVolume),
                ("", ThruSection), ("!", CreateShell),
                ("Protrude...", Extrude), ("", Revolve), ("!", Sweep),
                ("Fillet/Chamfer", Fillet), ("!", Chamfer),                                
                ("Copy/Remove...", Copy), ("", Remove), ("", Remove2), ("!", RemoveFaces),
                ("Translate...", Move,), ("", MoveByPoints), ("", Rotate), ("", RotateCenterPoints),
                ("", Flip), ("!", Scale),
                ("Array...", Array), ("", ArrayRot), ("", ArrayByPoints), ("", ArrayRotByPoints),
                ("!", ArrayPath), 
                ("Boolean...", Union), ("", MergeFace), ("", Intersection),
                ("", Difference), ("", Fragments), ("!", SplitByPlane),
                ("WorkPlane...", WorkPlane), ("", WorkPlaneByPoints),
                ("", WPParallelToPlane), ("!", WPNormalToPlane),
                ("Import...", BrepImport), ("", CADImport), ("!", healCAD),
                ]

    def get_special_menu(self, evt):
        if not hasattr(self, '_gso'):
            return [('Build All', self.onBuildAll, None),
                    ('Export Brep', self.onExportBrep, None)]

        dlg = evt.GetEventObject().GetTopLevelParent()
        viewer = dlg.GetParent()
        selection = viewer.dom_bdr_sel

        if np.sum([len(x) for x in selection]) == 0:
            menu = [('Build All', self.onBuildAll, None),
                    ('Export Brep', self.onExportBrep, None)]
        else:
            menu = [('Build All', self.onBuildAll, None),
                    ('Export Selected Entity',
                     self.onExportSelectedBrep, None),
                    ('Export Brep', self.onExportBrep, None)]

        if self._gso.child_alive():
            m2 = [('---', None, None),
                  ('Terminate geometry process',
                   self.onTerminateChild, None),]
            menu.extend(m2)
        return menu
