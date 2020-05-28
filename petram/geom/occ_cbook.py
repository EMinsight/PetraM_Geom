from petram.geom.geom_id import (GeomIDBase, VertexID, LineID, SurfaceID, VolumeID,
                                 LineLoopID, SurfaceLoopID)
import numpy as np

hasOCC = False

try:
    from OCC.Core.GeomAPI import GeomAPI_Interpolate
    from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
    from OCC.Core.TopLoc import TopLoc_Location
    from OCC.Core.TopExp import (TopExp_Explorer,
                                 topexp_MapShapes,
                                 topexp_MapShapesAndAncestors)
    from OCC.Core.BRep import BRep_Builder, BRep_Tool
    from OCC.Core.BRepTools import (breptools_Write,
                                    breptools_Read,
                                    breptools_Clean)
    from OCC.Core.TopTools import (TopTools_IndexedMapOfShape,
                                   TopTools_IndexedDataMapOfShapeListOfShape,
                                   TopTools_ListIteratorOfListOfShape,
                                   TopTools_ListOfShape)
    from OCC.Core.ShapeFix import (ShapeFix_Shape,
                                   ShapeFix_Solid,
                                   ShapeFix_Shell,
                                   ShapeFix_Face,
                                   ShapeFix_Wire,
                                   ShapeFix_Wireframe,
                                   ShapeFix_FixSmallFace)
    from OCC.Core.TopoDS import (TopoDS_CompSolid,
                                 TopoDS_Compound,
                                 TopoDS_Shape,
                                 TopoDS_Solid,
                                 TopoDS_Shell,
                                 TopoDS_Face,
                                 TopoDS_Wire,
                                 TopoDS_Edge,
                                 TopoDS_Vertex,
                                 topods_Compound,
                                 topods_CompSolid,
                                 topods_Solid,
                                 topods_Shell,
                                 topods_Face,
                                 topods_Wire,
                                 topods_Edge,
                                 topods_Vertex)
    from OCC.Core.TopAbs import (TopAbs_COMPSOLID,
                                 TopAbs_COMPOUND,                                 
                                 TopAbs_SOLID,
                                 TopAbs_SHELL,
                                 TopAbs_FACE,
                                 TopAbs_WIRE,
                                 TopAbs_EDGE,
                                 TopAbs_VERTEX)
    from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakePrism,
                                      BRepPrimAPI_MakeRevol,
                                      BRepPrimAPI_MakeCone,
                                      BRepPrimAPI_MakeWedge,
                                      BRepPrimAPI_MakeSphere,
                                      BRepPrimAPI_MakeTorus,
                                      BRepPrimAPI_MakeCylinder,)
    from OCC.Core.BRepFilletAPI import (BRepFilletAPI_MakeFillet,
                                        BRepFilletAPI_MakeChamfer)
    from OCC.Core.BRepOffsetAPI import (BRepOffsetAPI_MakePipe,
                                        BRepOffsetAPI_MakeFilling,
                                        BRepOffsetAPI_ThruSections,
                                        BRepOffsetAPI_NormalProjection)
    from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_Sewing,
                                         BRepBuilderAPI_Copy,
                                         BRepBuilderAPI_Transform,
                                         BRepBuilderAPI_GTransform,
                                         BRepBuilderAPI_MakeSolid,
                                         BRepBuilderAPI_MakeShell,
                                         BRepBuilderAPI_MakeFace,
                                         BRepBuilderAPI_MakeWire,
                                         BRepBuilderAPI_MakeEdge,
                                         BRepBuilderAPI_MakeVertex)
    from OCC.Core.BRepAlgoAPI import (BRepAlgoAPI_Fuse,
                                      BRepAlgoAPI_Cut,
                                      BRepAlgoAPI_Common,
                                      BRepAlgoAPI_BuilderAlgo,
                                      BRepAlgoAPI_Defeaturing)
    from OCC.Core.gp import (gp_Ax1, gp_Ax2, gp_Pnt,
                             gp_Dir, gp_Pnt2d, gp_Trsf,
                             gp_Vec, gp_XYZ, gp_GTrsf, gp_Mat)
    from OCC.Core.GC import (GC_MakeArcOfCircle,
                             GC_MakeSegment,
                             GC_MakeCircle)
    from OCC.Core.BOPTools import BOPTools_AlgoTools3D

    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import (brepgprop_LinearProperties,
                                    brepgprop_SurfaceProperties)

    __ex1 = TopExp_Explorer()
    __ex2 = TopExp_Explorer()
    
    __expparam = {'compound': (TopAbs_COMPOUND, topods_Compound, ''),
                  'compsolid': (TopAbs_COMPSOLID, topods_CompSolid, ''),
                  'compsolid': (TopAbs_SOLID, topods_Solid, ''),                  
                  'shell': (TopAbs_SHELL, topods_Shell, 'solid'),
                  'face': (TopAbs_FACE, topods_Face, 'shell'),
                  'wire': (TopAbs_WIRE, topods_Wire, 'face'),
                  'edge': (TopAbs_EDGE, topods_Edge, 'wire'),
                  'vertex': (TopAbs_VERTEX, topods_Vertex, 'edge')}
    __topo_names = ('solid', 'shell', 'face', 'wire', 'edge', 'vertex')
    hasOCC = True

except ImportError:
    import traceback
    print(" ****** OCC module import error ******")
    traceback.print_exc()


def iter_shape(shape, shape_type='shell', exclude_parent=False):
    '''
    iterate over shape. this allows to write

    for subshape in iter_shape(shape, 'shell'):
        ...
    '''
    args = [__expparam[shape_type][0], ]
    cast = __expparam[shape_type][1]

    if exclude_parent:
        args.append(__expparam[shape_type][2])

    __ex1.Init(shape, *args)
    while __ex1.More():
        sub_shape = cast(__ex1.Current())
        yield sub_shape
        __ex1.Next()


def iterdouble_shape(shape_in, inner_type='shell'):
    outer_type = __expparam[inner_type][2]

    args1 = [__expparam[outer_type][0], ]
    args2 = [__expparam[inner_type][0], ]
    cast1 = __expparam[outer_type][1]
    cast2 = __expparam[inner_type][1]

    __ex1.Init(shape_in, *args1)

    while __ex1.More():
        p_shape = cast1(__ex1.Current())
        __ex2.Init(p_shape, *args2)

        while __ex2.More():
            shape = cast2(__ex2.Current())
            yield shape, p_shape
            __ex2.Next()
        __ex1.Next()


def do_rotate(shape, ax, an, txt=''):
    ax = [float(x) for x in ax]
    trans = gp_Trsf()
    axis_revolution = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(ax[0], ax[1], ax[2]))
    trans.SetRotation(axis_revolution, an)
    transformer = BRepBuilderAPI_Transform(trans)
    transformer.Perform(shape)
    if not transformer.IsDone():
        assert False, "can not rotate (WP "+txt+")"
    return transformer.ModifiedShape(shape)


def do_translate(shape, delta):
    delta = [float(x) for x in delta]
    trans = gp_Trsf()
    trans.SetTranslation(gp_Vec(delta[0], delta[1], delta[2]))
    transformer = BRepBuilderAPI_Transform(trans)
    transformer.Perform(shape)
    if not transformer.IsDone():
        assert False, "can not translate (WP) "
    return transformer.ModifiedShape(shape)

def get_topo_list():
    lists = {}
    for x in __topo_names:
        cls = 'topo_list_'+x
        lists[x] = globals()[cls]()
    return lists

def calc_wp_projection(c1, a1, a2):
    x1 = np.array([1., 0., 0.])

    ax = np.cross(x1, a1)
    an = np.arctan2(np.sqrt(np.sum(ax**2)), np.dot(a1, x1))

    if np.sum(ax**2) == 0.0:
        if an != 0.0:
            # if a1 is [0, 0, -1], rotate 180 deg
            ax = np.array([0, 1, 0])
            an = np.pi
        else:
            ax = x1
            an = 0.0
    if np.sum(ax**2) != 0.0 and an != 0.0:
        ax1 = ax
        an1 = an
    else:
        ax1 = x1
        an1 = 0.0

    from petram.geom.geom_utils import rotation_mat
    R = rotation_mat(ax1, an1)
    y2 = np.dot(R, np.array([0, 1, 0]))
    ax = a1
    aaa = np.cross(a1, y2)
    an = np.arctan2(np.dot(a2, aaa), np.dot(a2, y2))

    if np.sum(ax**2) == 0.0 and an != 0.0:
        # rotate 180 deg around a1
        ax2 = a1
        an2 = np.pi
    else:
        ax2 = ax
        an2 = an

    if c1[0] != 0.0 or c1[1] != 0.0 or c1[2] != 0.0:
        cxyz = c1
    else:
        cxyz = None

    return ax1, an1, ax2, an2, cxyz

'''
def prep_maps(shape, return_all=True):
    solidMap = TopTools_IndexedMapOfShape()
    faceMap = TopTools_IndexedMapOfShape()
    edgeMap = TopTools_IndexedMapOfShape()
    vertMap = TopTools_IndexedMapOfShape()

    topexp_MapShapes(shape, TopAbs_SOLID, solidMap)
    topexp_MapShapes(shape, TopAbs_FACE, faceMap)
    topexp_MapShapes(shape, TopAbs_EDGE, edgeMap)
    topexp_MapShapes(shape, TopAbs_VERTEX, vertMap)

    if not return_all:
        return (solidMap, faceMap, edgeMap, vertMap)

    shellMap = TopTools_IndexedMapOfShape()
    wireMap = TopTools_IndexedMapOfShape()
    topexp_MapShapes(shape, TopAbs_SHELL, shellMap)
    topexp_MapShapes(shape, TopAbs_WIRE, wireMap)

    return (solidMap, shellMap, faceMap, wireMap,
            edgeMap, vertMap)
'''

def prep_maps(shape, return_all=True, return_compound=False):

    solidMap = TopTools_IndexedMapOfShape()
    faceMap = TopTools_IndexedMapOfShape()
    edgeMap = TopTools_IndexedMapOfShape()
    vertMap = TopTools_IndexedMapOfShape()

    topexp_MapShapes(shape, TopAbs_SOLID, solidMap)
    topexp_MapShapes(shape, TopAbs_FACE, faceMap)
    topexp_MapShapes(shape, TopAbs_EDGE, edgeMap)
    topexp_MapShapes(shape, TopAbs_VERTEX, vertMap)

    maps = {'solid': solidMap, 'face': faceMap, 'edge': edgeMap,
            'vertex': vertMap}

    if not return_all:
        return maps

    shellMap = TopTools_IndexedMapOfShape()
    wireMap = TopTools_IndexedMapOfShape()
    topexp_MapShapes(shape, TopAbs_SHELL, shellMap)
    topexp_MapShapes(shape, TopAbs_WIRE, wireMap)

    maps['shell'] = shellMap
    maps['wire'] = wireMap
    
    if not return_compound:
        return maps

    compoundMap = TopTools_IndexedMapOfShape()
    compsolidMap = TopTools_IndexedMapOfShape()
    topexp_MapShapes(shape, TopAbs_COMPOUND, compoundMap)
    topexp_MapShapes(shape, TopAbs_COMPSOLID, compsolidMap)

    maps['compound'] = compoundMap
    maps['compsolid'] = compsolidMap

    return maps

def register_shape(shape, topolists):
    maps = prep_maps(shape)
    seens = {x: topo_seen(mapping=maps[x]) for x in maps}

    new_objs = []
    # registor solid
    for solid in iterd_shape(shape, 'solid'):
        if seens['solid'].check_shape(solid) == 0:
            solid_id = topolists['solid'].add(solid)
            new_objs.append(solid_id)

    def register_topo(shape, shape_type, add_newobj=False):
        seen = seens[shape_type]
        topo_ll = topolists[shape_type]

        for sub_shape, sub_shape_p in iterdouble_shape(shape, shape_type):
            if seen.check_shape(sub_shape) == 0:
                topo_id = topo_ll.add(sub_shape)
        for sub_shape in iter_shape(shape, shape_type, exclude_parent=True):
            if seen.check_shape(sub_shape) == 0:
                topo_id = topo_ll.add(sub_shape)
                if add_newobj:
                    new_objs.append(topo_id)

    register_topo(shape, 'shell')
    register_topo(shape, 'face', True)
    register_topo(shape, 'wire')
    register_topo(shape, 'edge', True)
    register_topo(shape, 'vertex', True)

    '''
    ex1 = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex1.More():
        solid = topods_Solid(ex1.Current())
        if usolids.check_shape(solid) == 0:
             solid_id = self.solids.add(solid)
             new_objs.append(solid_id)
        ex1.Next()

    def register_topo(shape, ucounter, topabs, topabs_p, topods, topods_p,
                      topo_list, dim=-1):
        ex1 = TopExp_Explorer(shape, topabs_p)
        while ex1.More():
            topo_p = topods_p(ex1.Current())
            ex2 = TopExp_Explorer(topo_p, topabs)
            while ex2.More():
                topo = topods(ex2.Current())
                if ucounter.check_shape(topo) == 0:
                    topo_id = topo_list.add(topo)
                ex2.Next()
            ex1.Next()
        ex1.Init(shape, topabs, topabs_p)
        while ex1.More():
            topo = topods(ex1.Current())
            if ucounter.check_shape(topo) == 0:
                topo_id = topo_list.add(topo)
                if dim != -1:
                    new_objs.append(topo_id)
            ex1.Next()
    
    register_topo(shape, ushells, TopAbs_SHELL, TopAbs_SOLID,
                  topods_Shell, topods_Solid, self.shells)
    register_topo(shape, ufaces, TopAbs_FACE, TopAbs_SHELL,
                  topods_Face, topods_Shell, self.faces, dim=2)
    register_topo(shape, uwires, TopAbs_WIRE, TopAbs_FACE,
                  topods_Wire, topods_Face, self.wires)
    register_topo(shape, uedges, TopAbs_EDGE, TopAbs_WIRE,
                  topods_Edge, topods_Wire, self.edges, dim=1)
    register_topo(shape, uvertices, TopAbs_VERTEX, TopAbs_EDGE,
                  topods_Vertex, topods_Edge,self.vertices, dim=0)
    '''


class topo_seen(list):
    def __init__(self, mapping):
        self.mapping = mapping
        self.check = np.array([0] * mapping.Size())

    def check_shape(self, x):
        i = self.mapping.FindIndex(x) - 1
        ret = self.check[i]
        self.check[i] += 1
        return ret

    def seen(self, x):
        return self.check_shape(x) != 0

    def not_seen(self, x):
        return self.check_shape(x) == 0


class topo_list():
    name = 'base'
    myclass = type(None)

    def __init__(self):
        self.gg = {0: {}, }
        self.d = self.gg[0]
        self.next_id = 0

    def add(self, shape):
        if not isinstance(shape, self.myclass):
            assert False, ("invalid object type" + self.myclass.__name__ +
                           ':' + shape.__class__.__name__)
        self.next_id += 1
        self.d[self.next_id] = shape
        return self.next_id

    def new_group(self):
        group = max(list(self.gg.keys()))+1
        self.set_group(group)
        return group

    def set_group(self, group):
        if not group in self.gg:
            self.gg[group] = {}
        self.d = self.gg[group]

    def current_group(self):
        for g in self.gg:
            if g is self.d:
                return g

    def __iter__(self):
        return self.d.__iter__()

    def __len__(self):
        return len(self.d)

    def __getitem__(self, val):
        return self.d[int(val)]

    def __contains__(self, val):
        return val in self.d

    def get_item_from_group(self, val, group=0):
        return self.gg[group][val]

    def is_toplevel(self, *args):
        del args  # not used
        assert False, "subclass need to add this"

    def synchronize(self, mapper, action='remove', verbose=False):
        if verbose:
            print("Synchronize:", self.name, mapper.Size())

        if action in ('remove', 'both'):
            removal = []
            found_idx = []
            for k in self.d:
                shape = self.d[k]
                if not mapper.Contains(shape):
                    removal.append(k)
            for k in removal:
                del self.d[k]
            if verbose:
                print("removed gid", removal, list(self.d))

        if action in ('add', 'both'):
            found_idx = []
            new_gids = []
            for k in self.d:
                shape = self.d[k]
                if mapper.Contains(shape):
                    idx = mapper.FindIndex(shape)
                    found_idx.append(idx)
            tmp = np.arange(1, mapper.Size() + 1)
            new_shape_idx = tmp[np.in1d(tmp, np.array(found_idx), invert=True)]

            for idx in new_shape_idx:
                shape = mapper(int(idx))
                new_gids.append(self.add(shape))

            if verbose:
                print("added gid", new_gids)


class topo_list_vertex(topo_list):
    name = 'vertex'
    myclass = TopoDS_Vertex

    def child_generator(self, val):
        del val  # unused
        return []

    def is_toplevel(self, val, compound):
        mapper = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            compound, TopAbs_VERTEX, TopAbs_EDGE, mapper)
        shape = self[val]
        if mapper.FindFromKey(shape).Size() == 0:
            return True
        return False

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_VERTEX, mapper)
        return mapper

    def get_child_mapper(self, args):
        del args  # unused
        return None

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return VertexID(ret)

    def keys(self):
        return [VertexID(x) for x in self.d]


class topo_list_edge(topo_list):
    name = 'edge'
    myclass = TopoDS_Edge

    def get_children(self, val):
        shape = self[val]
        return iter_shape(shape, 'vertex')
        '''
        ex1 = TopExp_Explorer(shape, TopAbs_VERTEX)
        while ex1.More():
            vertex = topods_Vertex(ex1.Current())
            yield vertex
            ex1.Next()
        '''

    def is_toplevel(self, val, compound):
        mapper = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            compound, TopAbs_EDGE, TopAbs_FACE, mapper)
        shape = self[val]
        if mapper.FindFromKey(shape).Size() == 0:
            return True
        return False

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_EDGE, mapper)
        return mapper

    def get_chilld_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_VERTEX, mapper)
        return mapper

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return LineID(ret)

    def keys(self):
        return [LineID(x) for x in self.d]


class topo_list_wire(topo_list):
    name = 'wire'
    myclass = TopoDS_Wire

    def get_children(self, val):
        shape = self[val]
        return iter_shape(shape, 'edge')
        '''
        ex1 = TopExp_Explorer(shape, TopAbs_EDGE)
        while ex1.More():
            edge = topods_Edge(ex1.Current())
            yield edge
            ex1.Next()
        '''

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_WIRE, mapper)
        return mapper

    def get_chilld_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_EDGE, mapper)
        return mapper

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return LineLoopID(ret)

    def keys(self):
        return [LineLoopID(x) for x in self.d]


class topo_list_face(topo_list):
    name = 'face'
    myclass = TopoDS_Face

    def get_children(self, val):
        shape = self[val]
        return iter_shape(shape, 'wire')
        #ex1 = TopExp_Explorer(shape, TopAbs_WIRE)
        # while ex1.More():
        #    wire = topods_Wire(ex1.Current())
        #    yield wire
        #    ex1.Next()

    def is_toplevel(self, val, compound):
        mapper = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            compound, TopAbs_FACE, TopAbs_SOLID, mapper)
        shape = self[val]
        if mapper.FindFromKey(shape).Size() == 0:
            return True
        return False

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_FACE, mapper)
        return mapper

    def get_chilld_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_WIRE, mapper)
        return mapper

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return SurfaceID(ret)

    def keys(self):
        return [SurfaceID(x) for x in self.d]


class topo_list_shell(topo_list):
    name = 'shell'
    myclass = TopoDS_Shell

    def get_children(self, val):
        shape = self[val]
        return iter_shape(shape, 'face')
        '''
        ex1 = TopExp_Explorer(shape, TopAbs_FACE)
        while ex1.More():
            face = topods_Face(ex1.Current())
            yield face
            ex1.Next()
        '''

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_SHELL, mapper)
        return mapper

    def get_chilld_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_FACE, mapper)
        return mapper

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return SurfaceLoopID(ret)

    def keys(self):
        return [SurfaceLoopID(x) for x in self.d]


class topo_list_solid(topo_list):
    name = 'solid'
    myclass = TopoDS_Solid

    def get_children(self, val):
        shape = self[val]
        return iter_shape(shape, 'shell')
        '''
        ex1 = TopExp_Explorer(shape, TopAbs_SHELL)
        while ex1.More():
            shell = topods_Shell(ex1.Current())
            yield shell
            ex1.Next()
        '''

    def get_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_SOLID, mapper)
        return mapper

    def get_chilld_mapper(self, shape):
        mapper = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_SHELL, mapper)
        return mapper

    def add(self, shape):
        ret = topo_list.add(self, shape)
        return VolumeID(ret)

    def keys(self):
        return [VolumeID(x) for x in self.d]


