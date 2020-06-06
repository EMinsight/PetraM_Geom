from petram.geom.gmsh_geom_model import GmshGeom

class OCCGeom(GmshGeom):
    has_2nd_panel = False
    
    @classmethod        
    def fancy_menu_name(self):
        return 'OCC Geometry'
    @classmethod
    def fancy_tree_name(self):
        return 'OCCSequence'

    def attribute_set(self, v):
        v = super(OCCGeom, self).attribute_set(v)
        v['long_edge_thr'] = 0.35
        v['small_edge_thr'] = 0.01
        return v
        
    def panel1_param(self):
        import wx
        return [["", "Geometry model using OpenCascade", 2, None],
                #["PreviewAlgorith", "Automatic", 4, {"style": wx.CB_READONLY,
                #                                     "choices": ["Auto", "MeshAdpat",
                #                                                 "Delaunay", "Frontal"]}],
                ["Preview Resolution (linear)", self.small_edge_thr, 300, None],
                ["Preview Resolution (angle)", self.long_edge_thr, 300, None],                
                ["Preview #threads", self.maxthreads, 400, None],
                [None, self.occ_parallel, 3, {"text": "OCC parallel boolean"}],
                [None, self.skip_final_frag, 3, {
                    "text": "Skip fragmentationn"}],
                #[None, self.use_1d_preview, 3, {"text": "Use line preview"}],
                #[None, self.use_occ_preview, 3, {
                #    "text": "OCC preview (in dev.)"}],
                #[None, self.use_curvature, 3, {
                #    "text": "Consider curvature in preview generation"}],
                [None, None, 341, {"label": "Finalize Geom",
                                   "func": 'onBuildAll',
                                   "noexpand": True}], ]

    def get_panel1_value(self):
        return [None, self.small_edge_thr, self.long_edge_thr,
                self.maxthreads, self.occ_parallel, self.skip_final_frag, self]

    def import_panel1_value(self, v):
        self.small_edge_thr = float(v[1])        
        self.long_edge_thr = float(v[2])
        self.maxthreads = int(v[3])        
        self.occ_parallel = v[4]        
        self.skip_final_frag = v[5]
        self.use_occ_preview = True

    def get_possible_child(self):
        from petram.geom.gmsh_primitives import (Point, PointCenter, PointByUV, PointOnEdge,
                                                 PointCircleCenter, Line, Spline,
                                                 Circle, CircleByAxisPoint, CircleBy3Points,
                                                 Rect, Polygon, OCCPolygon, Box, Ball,
                                                 Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep,
                                                 LineLoop, CreateLine, CreateSurface, CreateVolume,
                                                 SurfaceLoop, Union, Union2, MergeFace, Intersection, Difference, Fragments,
                                                 SplitByPlane, Copy, Remove, Remove, Remove2, RemoveFaces,
                                                 Move, Rotate, Flip, Scale, WorkPlane,
                                                 WorkPlaneByPoints, healCAD, CADImport, BrepImport,
                                                 Fillet, Chamfer,
                                                 Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                                                 ThruSection, RotateCenterPoints, MoveByPoints, ExtendedLine)
        return [Point, PointCenter, PointOnEdge, PointByUV, PointCircleCenter,
                Line, Circle, CircleByAxisPoint, CircleBy3Points,
                Rect, Polygon,  OCCPolygon, Spline, Box,
                Ball, Cone, Wedge, Cylinder,
                Torus, CreateLine, CreateSurface, CreateVolume, LineLoop, SurfaceLoop,
                Extrude, Revolve, Sweep, Union, Union2, MergeFace,
                Intersection, Difference, Fragments, SplitByPlane, Copy, Remove,
                Remove2, RemoveFaces, Move, Rotate,
                Flip, Scale, WorkPlane, WorkPlaneByPoints, healCAD, CADImport, BrepImport,
                Fillet, Chamfer, Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                ThruSection, RotateCenterPoints, MoveByPoints, ExtendedLine]

    def get_possible_child_menu(self):
        from petram.geom.gmsh_primitives import (Point, PointCenter, PointCircleCenter,
                                                 PointOnEdge, PointByUV,  Line, Spline,  
                                                 Circle, CircleByAxisPoint, CircleBy3Points,
                                                 Rect, OCCPolygon, Box, Ball,
                                                 Cone, Wedge, Cylinder, Torus, Extrude, Revolve, Sweep,
                                                 LineLoop, CreateLine, CreateSurface, CreateVolume,
                                                 SurfaceLoop, Union, Union2, MergeFace, Intersection, Difference, Fragments,
                                                 SplitByPlane, Copy, Remove, Remove2, RemoveFaces,
                                                 Move, Rotate, Flip, Scale,
                                                 WorkPlane, WorkPlaneByPoints, healCAD, CADImport, BrepImport,
                                                 Fillet, Chamfer,
                                                 Array, ArrayRot, ArrayByPoints, ArrayRotByPoints,
                                                 ThruSection, RotateCenterPoints, MoveByPoints, ExtendedLine)
        return [("Points...", Point), ("", PointCenter), ("", PointOnEdge),
                ("", PointCircleCenter), ("!", PointByUV),
                ("Lines", Line), ("!", ExtendedLine),
                ("Circle...", Circle), ("", CircleByAxisPoint), ("!", CircleBy3Points),
                ("", Rect),
                ("", Spline), ("", Fillet), ("", Chamfer),
                ("3D shape...", Box),
                ("", Ball), ("", Cone), ("", Wedge), ("", Cylinder),
                ("!", Torus),
                ("Create...", CreateLine), ("", CreateSurface), ("", CreateVolume),
                ("", OCCPolygon), ("!", ThruSection), #("", SurfaceLoop),
                ("Protrude...", Extrude), ("", Revolve), ("!", Sweep),
                ("Copy/Remove...", Copy), ("", Remove), ("", Remove2), ("!", RemoveFaces),
                ("Translate...", Move,), ("", MoveByPoints), ("", Rotate), ("", RotateCenterPoints),
                ("", Flip), ("!", Scale),
                ("Array...", Array), ("", ArrayRot), ("", ArrayByPoints), ("!", ArrayRotByPoints), 
                ("Boolean...", Union), ("", MergeFace), ("", Intersection),
                ("", Difference), ("", Fragments), ("!", SplitByPlane),
                ("WorkPlane...", WorkPlane), ("!", WorkPlaneByPoints),
                ("Import...", BrepImport), ("!", CADImport)
                ]
        
