from __future__ import print_function
from petram.geom.geom_id import (GeomIDBase, VertexID, LineID, SurfaceID, VolumeID,
                                 LineLoopID, SurfaceLoopID)

from OCC.Core.GeomAPI import GeomAPI_Interpolate

from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopExp import (TopExp_Explorer,
                             topexp_MapShapes,
                             topexp_MapShapesAndAncestors)
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.BRepTools import breptools_Write
from OCC.Core.TopTools import (TopTools_IndexedMapOfShape,
                               TopTools_IndexedDataMapOfShapeListOfShape,
                               TopTools_ListIteratorOfListOfShape,
                               TopTools_ListOfShape)
from OCC.Core.ShapeFix import (ShapeFix_Solid,
                               ShapeFix_Shell,
                               ShapeFix_Face)
from OCC.Core.TopoDS import (TopoDS_Compound,
                             TopoDS_Shape,
                             TopoDS_Solid,
                             TopoDS_Shell,
                             TopoDS_Face,
                             TopoDS_Wire,
                             TopoDS_Edge,
                             TopoDS_Vertex,
                             topods_Solid,
                             topods_Shell,
                             topods_Face,
                             topods_Wire,
                             topods_Edge,
                             topods_Vertex)
from OCC.Core.TopAbs import (TopAbs_SOLID,
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
                                    BRepOffsetAPI_MakeFilling)
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
from OCC.Core.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.Core.BOPTools import BOPTools_AlgoTools3D

from petram.geom.gmsh_geom_model import get_geom_key
from petram.phys.vtable import VtableElement, Vtable
import gmsh

import os
import numpy as np
import time
import tempfile
from collections import defaultdict
import multiprocessing as mp
from six.moves.queue import Empty as QueueEmpty

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('OCCGeomWrapper')


class Counter():
    def __init__(self):
        self.value = 0

    def increment(self, x):
        self.value = self.value + x

    def __call__(self):
        return self.value


def id2dimtag(en):
    if isinstance(en, VertexID):
        return (0, int(en))

    if isinstance(en, LineID):
        return (1, int(en))

    if isinstance(en, SurfaceID):
        return (2, int(en))

    if hasattr(en, 'surface'):
        return (2, int(en))

    if isinstance(en, VolumeID):
        return (3, int(en))

    assert False, "Illegal entity" + str(en)
    #    return None


def get_dimtag(entity):
    dimtags = []
    for en in entity:
        dimtags.append(id2dimtag(en))
    return dimtags


def dimtag2id(dimtags):
    out3 = []
    out2 = []
    out1 = []
    out0 = []
    for dim, tag in dimtags:
        if dim == 3 and not tag in out3:
            out3.append(VolumeID(tag))
        elif dim == 2 and not tag in out2:
            out2.append(SurfaceID(tag))
        elif dim == 1 and not tag in out1:
            out1.append(LineID(tag))
        elif dim == 0 and not tag in out1:
            out0.append(VertexID(tag))
    return out3 + out2 + out1 + out0


def find_combined_bbox(model, dimtags):
    xmax = -np.inf
    xmin = np.inf
    ymax = -np.inf
    ymin = np.inf
    zmax = -np.inf
    zmin = np.inf

    def update_maxmin(dim, tag, xmin, ymin, zmin, xmax, ymax, zmax):
        x1, y1, z1, x2, y2, z2 = model.getBoundingBox(dim, tag)
        xmax = np.max([xmax, x2])
        ymax = np.max([ymax, y2])
        zmax = np.max([zmax, z2])
        xmin = np.min([xmin, x1])
        ymin = np.min([ymin, y1])
        zmin = np.min([zmin, z1])
        return xmin, ymin, zmin, xmax, ymax, zmax

    for dim, tag in dimtags:
        xmin, ymin, zmin, xmax, ymax, zmax = update_maxmin(dim, tag,
                                                           xmin, ymin, zmin,
                                                           xmax, ymax, zmax)
    return xmin, ymin, zmin, xmax, ymax, zmax


class topo_seen(list):
    def __init__(self, mapping):
        self.mapping = mapping
        self.check = np.array([0] * mapping.Size())

    def check_shape(self, x):
        i = self.mapping.FindIndex(x) - 1
        ret = self.check[i]
        self.check[i] += 1
        return ret


class topo2id():
    def __init__(self, dd, mapper):
        self.mapper = mapper
        self.mapperout2k = {mapper.FindIndex(dd[k]): k for k in dd}
        #self._d = [(dd[k], k) for k in dd]

    def __getitem__(self, val):
        out = self.mapper.FindIndex(val)
        return self.mapperout2k[out]
        #assert False, "ID is not found in ap from Topo to ID"


class topo_list():
    name = 'base'
    myclass = type(None)
    def __init__(self):
        self.d = {}
        self.next_id = 0

    def add(self, shape):
        if not isinstance(shape, self.myclass):
            assert False, ("invalid object type" + self.myclass.__name__ +
                           ':' + shape.__class__.__name__)
        self.next_id += 1
        self.d[self.next_id] = shape
        return self.next_id

    def __iter__(self):
        return self.d.__iter__()

    def __len__(self):
        return len(self.d)

    def __getitem__(self, val):
        return self.d[int(val)]

    def __contains__(self, val):
        return val in self.d

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
            if verbose: print("removed gid", removal)

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
        ex1 = TopExp_Explorer(shape, TopAbs_VERTEX)
        while ex1.More():
            vertex = topods_Vertex(ex1.Current())
            yield vertex
            ex1.Next()

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
        ex1 = TopExp_Explorer(shape, TopAbs_EDGE)
        while ex1.More():
            edge = topods_Edge(ex1.Current())
            yield edge
            ex1.Next()

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
        ex1 = TopExp_Explorer(shape, TopAbs_WIRE)
        while ex1.More():
            wire = topods_Wire(ex1.Current())
            yield wire
            ex1.Next()

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
        ex1 = TopExp_Explorer(shape, TopAbs_FACE)
        while ex1.More():
            face = topods_Face(ex1.Current())
            yield face
            ex1.Next()

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
        ex1 = TopExp_Explorer(shape, TopAbs_SHELL)
        while ex1.More():
            shell = topods_Shell(ex1.Current())
            yield shell
            ex1.Next()

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


class Geometry():
    def __init__(self, **kwargs):
        self._point_loc = {}

        self.process_kwargs(kwargs)

        self.builder = BRep_Builder()
        self.bt = BRep_Tool()

        self.geom_sequence = []
        self._point = {}
        self._point_mask = []
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

        write_log = kwargs.pop('write_log', False)
        if write_log:
            self.logfile = tempfile.NamedTemporaryFile('w', delete=True)
        else:
            self.logfile = None

        self.queue = kwargs.pop("queue", None)
        self.p = None

    def process_kwargs(self, kwargs):
        self.geom_prev_res = kwargs.pop('PreviewResolution', 30)
        self.geom_prev_algorithm = kwargs.pop('PreviewAlgorithm', 2)
        self.occ_parallel = kwargs.pop('OCCParallel', 0)
        self.occ_boolean_tolerance = kwargs.pop('OCCBooleanTol', 1e-8)
        self.occ_geom_tolerance = kwargs.pop('OCCGeomTol', 1e-6)        
        self.maxthreads = kwargs.pop('Maxthreads', 1)
        self.skip_final_frag = kwargs.pop('SkipFrag', False)
        self.use_1d_preview = kwargs.pop('Use1DPreview', False)
        self.use_occ_preview = kwargs.pop('UseOCCPreview', False)
        self.long_edge_thr = kwargs.pop('LongEdgeThr', 0.1)
        self.small_edge_thr = kwargs.pop('SmallEdgeThr', 0.001)
        self.small_edge_seg = kwargs.pop('SmallEdgeSeg', 3)
        self.max_seg = kwargs.pop('MaxSeg', 30)

    def prep_topo_list(self):
        self.vertices = topo_list_vertex()
        self.edges = topo_list_edge()
        self.wires = topo_list_wire()
        self.faces = topo_list_face()
        self.shells = topo_list_shell()
        self.solids = topo_list_solid()
        
    def get_topo_list_for_gid(self, gid, child=0):
        ll = [self.vertices, self.edges, self.wires,
              self.faces, self.shells, self.solids]
        idx = gid.idx
        idx = idx - child
        if idx < 0:
            return None
        return ll[idx]
    
    def add_to_topo_list(self, shape):
        '''
        add shpae those kind is not known
        '''
        if isinstance(shape, TopoDS_Solid):
            gid = self.solids.add(shape)
        elif isinstance(shape, TopoDS_Shell):
            gid =self.shells.add(shape)
        elif isinstance(shape, TopoDS_Face):
            gid = self.faces.add(shape)
        elif isinstance(shape, TopoDS_Wire):
            gid = self.wires.add(shape)
        elif isinstance(shape, TopoDS_Edge):
            gid =self.edges.add(shape)
        elif isinstance(shape, TopoDS_Vertex):
            gid = self.vertices.add(shape)
        else:
            assert False, "Unkown shape type: " + type(shape)
        return gid
        
    def gid2shape(self, gid):
        d = self.get_topo_list_for_gid(gid)
        return d[int(gid)]

    def print_number_of_topo_objects(self):
        solidMap, faceMap, edgeMap, vertMap = self.prep_maps(
            self.shape, return_all=False)

        dprint1("Entity counts: solid/face/edge/vert : ",
                solidMap.Size(), faceMap.Size(), edgeMap.Size(), vertMap.Size())
        
    def count_topos(self):
        solidMap, faceMap, edgeMap, vertMap = self.prep_maps(
            self.shape, return_all=False)
        return (solidMap.Size(), faceMap.Size(), edgeMap.Size(), vertMap.Size())

    def bounding_box(self, shape=None, tolerance=1e-5):
        from OCC.Core.Bnd import Bnd_Box
        from OCC.Core.BRepBndLib import brepbndlib_Add

        shape = self.shape if shape is None else shape

        bbox = Bnd_Box()
        bbox.SetGap(tolerance)
        brepbndlib_Add(shape, bbox)
        values = bbox.Get()
        return [values[i] for i in range(6)]

    def get_esize(self):
        esize = {}
        for iedge in self.edges:
            x1, y1, z1, x2, y2, z2 = self.bounding_box(self.edges[iedge])
            s = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5

            esize[iedge] = s
        return esize

    def get_vcl(self, l, esize):
        lcar = defaultdict(lambda: np.inf)
        for iedge in esize:
            if not iedge in l:
                continue
            iverts = l[iedge]
            for ivert in iverts:
                lcar[ivert] = min(lcar[ivert], esize[iedge])

        return dict(lcar)

    def get_target1(self, objs, targets, cls):
        # this is when target type is given
        if cls == 'l':
            cc = LineID
        elif cls == 'v':
            cc = VolumeID
        elif cls == 'f':
            cc = SurfaceID
        elif cls == 'p':
            cc = VertexID
        else:
            pass
        gids = [objs[t] if t in objs else cc(t) for t in targets]
        return gids

    def get_target2(self, objs, targets):
        # this is when target type is not given
        ret = []

        for t in targets:
            if t in objs:
                ret.append(objs[t])
            else:
                if t.startswith("p"):
                    ret.append(VertexID(int(t[1:])))
                if t.startswith("l"):
                    ret.append(LineID(int(t[1:])))
                if t.startswith("f"):
                    ret.append(SurfaceID(int(t[1:])))
                if t.startswith("v"):
                    ret.append(VolumeID(int(t[1:])))

        if len(ret) == 0:
            assert False, "empty imput objects: "+','.join(targets)

        return ret

    def get_point_coord(self, gid):
        if not gid in self.vertices:
            assert False, "can not find point: "+str(int(gid))
        shape = self.vertices[gid]
        pnt = self.bt.Pnt(shape)
        return np.array((pnt.X(), pnt.Y(), pnt.Z(),))

    def get_planesurface_normal(self, gid):
        '''
        return normal vector of flat surface and a representative point on 
        the plane
        '''
        if not gid in self.faces:
            assert False, "can not find surface: "+str(int(gid))

        shape = self.faces[gid]
        surface = self.bt.Surface(shape)

        uMin, uMax, vMin, vMax = surface.Bounds()

        dirc = gp_Dir()
        tool = BOPTools_AlgoTools3D()
        tool.GetNormalToSurface(surface, uMin, vMin, dirc)
        n1 = (dirc.X(), dirc.Y(), dirc.Z())
        tool.GetNormalToSurface(surface, uMin, vMax, dirc)
        n2 = (dirc.X(), dirc.Y(), dirc.Z())
        tool.GetNormalToSurface(surface, uMax, vMin, dirc)
        n3 = (dirc.X(), dirc.Y(), dirc.Z())

        if n1 != n2:
            assert False, "surface is not flat"
        if n1 != n3:
            assert False, "surface is not flat"

        ptx = gp_Pnt()
        surface.D0(uMin, vMin, ptx)
        ptx = (ptx.X(), ptx.Y(), ptx.Z())

        return np.array(n1), np.array(ptx)

    def write_brep(self, filename):

        comp = TopoDS_Compound()
        b = self.builder
        b.MakeCompound(comp)
        ex1 = TopExp_Explorer(self.shape, TopAbs_SOLID)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(self.shape, TopAbs_FACE, TopAbs_SHELL)
        while ex1.More():
            print("write face")
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(self.shape, TopAbs_EDGE, TopAbs_WIRE)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(self.shape, TopAbs_VERTEX, TopAbs_EDGE)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()

        dprint1("exproted brep file:", filename)
        breptools_Write(comp, filename)

    def prep_maps(self, shape, return_all=True):
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

    def synchronize_topo_list(self, **kwargs):
        solidMap, shellMap, faceMap, wireMap, edgeMap, vertMap = self.prep_maps(
            self.shape)
        self.solids.synchronize(solidMap, **kwargs)
        self.shells.synchronize(shellMap, **kwargs)
        self.faces.synchronize(faceMap, **kwargs)
        self.wires.synchronize(wireMap, **kwargs)
        self.edges.synchronize(edgeMap, **kwargs)
        self.vertices.synchronize(vertMap, **kwargs)


    @property
    def dim(self):
        if len(self.model.getEntities(3)) > 0:
            return 3
        if len(self.model.getEntities(2)) > 0:
            return 2
        if len(self.model.getEntities(1)) > 0:
            return 1
        return 0

    def add_point(self, p, lcar=0.0, mask=True):
        p = BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])).Shape()
        return self.vertices.add(p)

    def add_line(self, p1, p2):
        edgeMaker = BRepBuilderAPI_MakeEdge(
            self.vertices[p1], self.vertices[p2])
        edgeMaker.Build()
        if not edgeMaker.IsDone():
            assert False, "Can not make line"
        edge = edgeMaker.Edge()
        return self.edges.add(edge)

    def add_circle_arc(self, p1, p3, p2):
        bt = self.bt
        pnt1 = bt.Pnt(self.vertices[p1])
        pnt2 = bt.Pnt(self.vertices[p2])
        pnt3 = bt.Pnt(self.vertices[p3])        

        arc = GC_MakeArcOfCircle(pnt1, pnt3, pnt2)
        
        edgeMaker = BRepBuilderAPI_MakeEdge(arc.Value())
        edgeMaker.Build()
        if not edgeMaker.IsDone():
            assert False, "Can not make circle arc"
        edge = edgeMaker.Edge()

        return self.edges.add(edge)        

    def add_ellipse_arc(self, startTag, centerTag, endTag):
        a = self._point[startTag] - self._point[centerTag]
        b = self._point[endTag] - self._point[centerTag]
        if np.sum(a * a) > np.sum(b * b):
            l = self.factory.addEllipseArc(startTag, centerTag, endTag)
        else:
            l = self.factory.addEllipseArc(endTag, centerTag, startTag)
        return LineID(l)

    def add_spline(self, pos, tolerance=1e-5, periodic=False):
        from OCC.Core.TColgp import TColgp_HArray1OfPnt

        bt = self.bt

        pts = [BRepBuilderAPI_MakeVertex(gp_Pnt(p[0], p[1], p[2])).Shape()
               for p in pos]

        array = TColgp_HArray1OfPnt(1, len(pts))
        for i, p in enumerate(pts):
            array.SetValue(i + 1, bt.Pnt(p))

        itp = GeomAPI_Interpolate(array, periodic, tolerance)
        itp.Perform()
        if not itp.IsDone():
            assert False, "Can not interpolate points (add_spline)"

        start = pts[0]
        end = pts[-1]
        if periodic:
            edgeMaker = BRepBuilderAPI_MakeEdge(itp.Curve(), start, start)
        else:
            edgeMaker = BRepBuilderAPI_MakeEdge(itp.Curve(), start, end)
        edgeMaker.Build()
        if not edgeMaker.IsDone():
            assert False, "Can not make spline"
        edge = edgeMaker.Edge()

        l_id = self.edges.add(edge)
        if periodic:
            self.vertices.add(start)
        else:
            self.vertices.add(start)
            self.vertices.add(end)
        return l_id

    def add_plane_surface(self, tag):

        wire = self.wires[tag]
        faceMaker = BRepBuilderAPI_MakeFace(wire)
        faceMaker.Build()

        if not faceMaker.IsDone():
            assert False, "can not create face"

        face = faceMaker.Face()

        fixer = ShapeFix_Face(face)
        fixer.Perform()
        face = fixer.Face()
        f_id = self.faces.add(face)

        return f_id

    def add_plate_surface(self, gids_edge, gids_vertex):
        from OCC.Core.GeomPlate import (GeomPlate_BuildPlateSurface,
                                        GeomPlate_PointConstraint,
   	                                GeomPlate_MakeApprox)
        from OCC.Core.BRepTools import BRepTools_WireExplorer        
        from OCC.Core.BRepAdaptor import BRepAdaptor_HCurve
        from OCC.Core.BRepFill import BRepFill_CurveConstraint
        from OCC.Core.ShapeFix import ShapeFix_Face


        bt = BRep_Tool()
        BPSurf = GeomPlate_BuildPlateSurface(2, 150, 10)

        # make wire first
        wireMaker = BRepBuilderAPI_MakeWire()
        for gid in gids_edge:
            edge = self.edges[gid]
            wireMaker.Add(edge)
        wireMaker.Build()

        if not wireMaker.IsDone():
            assert False, "Failed to make wire"
        wire = wireMaker.Wire()


        # make wire constraints
        ex1 = BRepTools_WireExplorer(wire)
        while ex1.More():
            edge = topods_Edge(ex1.Current())
            C = BRepAdaptor_HCurve()
            C.ChangeCurve().Initialize(edge)
            Cont = BRepFill_CurveConstraint(C, 0)
            BPSurf.Add(Cont)
            ex1.Next()

        # make point constraints
        for gid in gids_vertex:
            vertex = self.vertices[gid]
            Pcont = GeomPlate_PointConstraint(bt.Pnt(vertex), 0)
            BPSurf.Add(Pcont)

        BPSurf.Perform()

        MaxSeg = 9
        MaxDegree = 8
        CritOrder = 0

        PSurf = BPSurf.Surface()

        dmax = max(0.0001, 10 * BPSurf.G0Error())
        Tol = 0.0001

        Mapp = GeomPlate_MakeApprox(PSurf, Tol, MaxSeg, MaxDegree,
                                    dmax, CritOrder)
        Surf = Mapp.Surface()
        uMin, uMax, vMin, vMax = Surf.Bounds()

        faceMaker = BRepBuilderAPI_MakeFace(Surf, uMin, uMax, vMin, vMax, 1e-6)
        result = faceMaker.Face()
        
        fix = ShapeFix_Face(result)
        fix.SetPrecision(self.occ_geom_tolerance)
        fix.Perform()
        fix.FixOrientation()
        result = fix.Face()

        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs

    def add_surface_filling(self, gids_edge, gids_vertex):
        from OCC.Core.GeomAbs import GeomAbs_C0
        from OCC.Core.BRepTools import BRepTools_WireExplorer

        bt = BRep_Tool()
        f = BRepOffsetAPI_MakeFilling()

        # make wire first
        wireMaker = BRepBuilderAPI_MakeWire()
        for gid in gids_edge:
            edge = self.edges[gid]
            wireMaker.Add(edge)
        wireMaker.Build()

        if not wireMaker.IsDone():
            assert False, "Failed to make wire"
        wire = wireMaker.Wire()
 
        #make wire constraints
        ex1 = BRepTools_WireExplorer(wire)
        while ex1.More():
            edge = topods_Edge(ex1.Current())
            f.Add(edge, GeomAbs_C0)
            ex1.Next()

        for gid in gids_vertex:
            vertex = self.vertices[gid]
            pnt = bt.Pnt(vertex)
            f.Add(pnt)

        f.Build()

        if not f.IsDone():
            assert False, "Cannot make filling"

        face = f.Shape()
        s = bt.Surface(face)

        faceMaker = BRepBuilderAPI_MakeFace(s, wire)
        self.wires.add(wire)

        result = faceMaker.Face()
        
        fix = ShapeFix_Face(result)
        fix.SetPrecision(self.occ_geom_tolerance)
        fix.Perform()
        fix.FixOrientation()
        result = fix.Face()

        face_id = self.faces.add(result)

        return face_id
    
    def add_line_loop(self, pts, sign=None):
        tags = list(np.atleast_1d(pts))

        wireMaker = BRepBuilderAPI_MakeWire()
        for t in tags:
            edge = self.edges[t]
            wireMaker.Add(edge)
        wireMaker.Build()

        if not wireMaker.IsDone():
            assert False, "Failed to make wire"
        wire = wireMaker.Wire()

        w_id = self.wires.add(wire)
        return w_id

    def add_curve_loop(self, pts):
        return self.add_line_loop(pts)

    def add_surface_loop(self, sl):
        tags = list(np.atleast_1d(sl))

        try:
            sewingMaker = BRepBuilderAPI_Sewing()
            for t in tags:
                face = self.faces[t]
                sewingMaker.Add(face)
            sewingMaker.Perform()
            result = sewingMaker.SewedShape()
        except BaseException:
            assert False, "Failed to sew faces"

        ex1 = TopExp_Explorer(result, TopAbs_SHELL)
        while ex1.More():
            shell = topods_Shell(ex1.Current())
            fixer = ShapeFix_Shell(shell)
            fixer.Perform()
            shell = fixer.Shell()
            break
            ex.Next()

        shell_id = self.shells.add(shell)

        return shell_id

    def add_volume(self, shells):
        tags = list(np.atleast_1d(shells))

        solidMaker = BRepBuilderAPI_MakeSolid()
        for t in tags:
            shell = self.shells[t]
            solidMaker.Add(shell)
        result = solidMaker.Solid()

        if not solidMaker.IsDone():
            assert False, "Failed to make solid"

        fixer = ShapeFix_Solid(result)
        fixer.Perform()
        result = topods_Solid(fixer.Solid())

        solid_id = self.solids.add(result)
        return solid_id

    def add_sphere(self, xyzc, radius, angle1, angle2, angle3):
        if radius <= 0:
            assert False, "Sphere radius should be > 0"

        if (angle3 <= 0 or angle3 > 2 * np.pi):
            assert False, "Cannot build sphere with angle <= 0 or angle > 2*pi"

        pnt = gp_Pnt(xyzc[0], xyzc[1], xyzc[2])
        s = BRepPrimAPI_MakeSphere(pnt, radius, angle1, angle2, angle3)

        s.Build()
        if not s.IsDone():
            assert False, "Could not create sphere"

        result = topods_Solid(s.Shape())
        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)
        
        return new_objs

    def add_cone(self, x, y, z, dx, dy, dz, r1, r2, angle):
        H = np.sqrt(dx * dx + dy * dy + dz * dz);
        if H == 0:
            assert False, "Cone hight must be > 0"
        if angle <= 0 :
            assert False, "Cone angle should be positive"

        pnt = gp_Pnt(x, y, z)
        vec = gp_Dir(dx/H, dy/H, dz/H)
        axis = gp_Ax2(pnt, vec)

        c = BRepPrimAPI_MakeCone(axis, r1, r2, H, angle)
        c.Build()
        if not c.IsDone():
            assert False, "Could not create cone"

        result = topods_Solid(c.Shape())
        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)
        
        return new_objs
        
    def add_wedge(self, xyz, dxyz, ltx):
        x, y, z = xyz
        dx, dy, dz = dxyz
        pnt = gp_Pnt(x, y, z)
        vec = gp_Dir(0, 0, 1)
        axis = gp_Ax2(pnt, vec)

        w = BRepPrimAPI_MakeWedge(axis, dx, dy, dz, ltx)
        w.Build()
        if not w.IsDone():
            assert False, "Could not create wedge"

        result = topods_Solid(w.Shape())
        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)
        
        return new_objs

    def add_cylinder(self, xyz, dxyz, r, angle):
        x, y, z = xyz
        dx, dy, dz = dxyz
        H = np.sqrt(dx * dx + dy * dy + dz * dz)
        if H == 0:
            assert False, "Cylinder height must be > 0"
        if r <= 0:
            assert False, "Cylinder radius must be > 0"
        if (angle <= 0 or angle > 2 * np.pi):
            assert False, "Cannot build a cylinder with angle <= 0 or angle > 2*pi"

        pnt = gp_Pnt(x, y, z)
        vec = gp_Dir(dx/H, dy/H, dz/H)
        axis = gp_Ax2(pnt, vec)

        cl = BRepPrimAPI_MakeCylinder(axis, r, H, angle)
        cl.Build()
        if not cl.IsDone():
            assert False, "Can not create cylinder"

        result = topods_Solid(cl.Shape())
        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)
        
        return new_objs

    def add_torus(self, xyz, r1, r2, angle):
        x, y, z = xyz
        pnt = gp_Pnt(x, y, z)
        vec = gp_Dir(0, 0, 1)
        axis = gp_Ax2(pnt, vec)

        if r1 <= 0:
            assert False, "Torus major radius must be > 0"
        if r2 <= 0:
            assert False, "Torus minor radius must be > 0"
        if (angle <= 0 or angle > 2 * np.pi):
            assert False, "Torus angle must be between 0, and 2*pi"

        t = BRepPrimAPI_MakeTorus(axis, r1, r2, angle)
        t.Build()
        if not t.IsDone():
            assert False, "Could not create torus"

        result = topods_Solid(t.Shape())
        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs

    def fillet(self, gid_vols, gid_curves, radii):

        comp = TopoDS_Compound()
        self.builder.MakeCompound(comp)


        for gid in gid_vols:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            self.builder.Add(comp, shape)
            self.remove(gid, recursive=True)            

        f = BRepFilletAPI_MakeFillet(comp)

        for kk, gid in enumerate(gid_curves):
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]

            if len(radii) == 1:
                f.Add(radii[0], shape)
            elif len(gid_curves) == len(radii):
                f.Add(radii[kk], shape)
            elif len(gid_curves)+1 == len(radii):
                f.Add(radii[kk], radii[kk+1], shape)
            else:
                assert False, "Wrong radius setting"

        f.Build()
        if not f.IsDone():
            assert False, "Can not make fillet"

        result = f.Shape()

        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs

    def chamfer(self, gid_vols, gid_curves, gid_faces, distances):

        comp = TopoDS_Compound()
        self.builder.MakeCompound(comp)


        for gid in gid_vols:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            self.builder.Add(comp, shape)
            self.remove(gid, recursive=True)
            
        f = BRepFilletAPI_MakeChamfer(comp)

        kk = 0
        for gid_c, gid_f in zip(gid_curves, gid_faces):
            topolist = self.get_topo_list_for_gid(gid_c)
            edge = topolist[gid_c]
            topolist = self.get_topo_list_for_gid(gid_f)
            face = topolist[gid_f]

            if len(distances) == 1:
                f.Add(distances[0], distances[0],  edge, face)
            elif len(distances) == len(gid_curves):
                f.Add(distances[kk], edge, face)
            elif len(distances) == len(gid_curves)*2:
                f.Add(distances[2*kk], distances[2*kk+1], edge, face)
            kk = kk+1

        f.Build()
        if not f.IsDone():
            assert False, "Can not make chamfer"

        result = f.Shape()

        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs


    def add_box(self, points):
        p1, p2, p3, p4, p5, p6, p7, p8 = points
        lcar = 0.0

        p1 = self.add_point(p1)
        p2 = self.add_point(p2)
        p3 = self.add_point(p3)
        p4 = self.add_point(p4)
        p5 = self.add_point(p5)
        p6 = self.add_point(p6)
        p7 = self.add_point(p7)
        p8 = self.add_point(p8)

        l1 = self.add_line(p1, p2)
        l2 = self.add_line(p2, p5)
        l3 = self.add_line(p5, p3)
        l4 = self.add_line(p3, p1)
        l5 = self.add_line(p1, p4)
        l6 = self.add_line(p2, p7)
        l7 = self.add_line(p5, p8)
        l8 = self.add_line(p3, p6)
        l9 = self.add_line(p4, p7)
        l10 = self.add_line(p7, p8)
        l11 = self.add_line(p8, p6)
        l12 = self.add_line(p6, p4)

        ll1 = self.add_curve_loop([l1, l2, l3, l4])
        ll2 = self.add_curve_loop([l5, l9, l6, l1])
        ll3 = self.add_curve_loop([l6, l10, l7, l2])
        ll4 = self.add_curve_loop([l7, l11, l8, l3])
        ll5 = self.add_curve_loop([l8, l12, l5, l4])
        ll6 = self.add_curve_loop([l9, l10, l11, l12])

        rec1 = self.add_plane_surface(ll1)
        rec2 = self.add_plane_surface(ll2)
        rec3 = self.add_plane_surface(ll3)
        rec4 = self.add_plane_surface(ll4)
        rec5 = self.add_plane_surface(ll5)
        rec6 = self.add_plane_surface(ll6)

        sl = self.add_surface_loop([rec1, rec2, rec3, rec4, rec5, rec6])

        v1 = self.add_volume(sl)

        return v1

    def update_topo_list_from_history(self, operator, list_of_shapes):
        iterator = TopTools_ListIteratorOfListOfShape(list_of_shapes)
        while iterator.More():
            shape = iterator.Value()            
            iterator.Next()
            # Do I need to do something with modified?
            # operator.Modified(shape)

            shape_gone = operator.IsDeleted(shape)
            if shape_gone:
                print("shape gone", shape_gone)

            shapes_new = operator.Generated(shape)
            iterator2 = TopTools_ListIteratorOfListOfShape(shapes_new)            
            while iterator2.More():
                shape_new = iterator2.Value()
                iterator2.Next()
                print("shape new", shape_new)

    def do_boolean(self, operation, gid_objs, gid_tools,
                   remove_tool=True, remove_obj=True,
                   keep_highest=False):

        if operation == 'fuse':
            operator = BRepAlgoAPI_Fuse()
        elif operation == 'cut':
            operator = BRepAlgoAPI_Cut()
        elif operation == 'common':
            operator = BRepAlgoAPI_Common()
        elif operation == 'fragments':
            operator = BRepAlgoAPI_BuilderAlgo()
        else:
            assert False, "Unknown boolean operation"

        objs = TopTools_ListOfShape()
        tools = TopTools_ListOfShape()

        for gid in gid_tools:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            tools.Append(shape)

        for gid in gid_objs:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            objs.Append(shape)


        operator.SetRunParallel(self.occ_parallel)
        operator.SetArguments(objs)

        if operation == 'fragments':
            tools.Clear()
        else:
            operator.SetTools(tools)

        if self.occ_boolean_tolerance > 0:
            operator.SetFuzzyValue(self.occ_boolean_tolerance)

        operator.Build()
        if not operator.IsDone():
            assert False, "boolean operation failed:" + operation

        result = operator.Shape()

        self.update_topo_list_from_history(operator, objs)
        self.update_topo_list_from_history(operator, tools)

        if remove_tool:
            for gid in gid_tools:
                self.remove(gid)
        if remove_obj:
            for gid in gid_objs:
                self.remove(gid)

        if keep_highest:
            result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs

    def union(self, gid_objs, gid_tools, remove_tool=True, remove_obj=True,
              keep_highest=False):

        return self.do_boolean('fuse', gid_objs, gid_tools,
                               remove_tool=remove_tool,
                               remove_obj=remove_obj,
                               keep_highest=keep_highest)
                               
    def intersection(self, gid_objs, gid_tools, remove_tool=True, remove_obj=True,
                     keep_highest=False):

        return self.do_boolean('common', gid_objs, gid_tools,
                               remove_tool=remove_tool,
                               remove_obj=remove_obj,
                               keep_highest=keep_highest)
                               
    def difference(self, gid_objs, gid_tools, remove_tool=True, remove_obj=True,
                   keep_highest=False):        
        
        return self.do_boolean('cut', gid_objs, gid_tools,
                               remove_tool=remove_tool,
                               remove_obj=remove_obj,
                               keep_highest=keep_highest)

    def fragments(self, gid_objs, gid_tools, remove_tool=True, remove_obj=True,
                 keep_highest=False):           

        gid_objs = gid_objs + gid_tools
        return self.do_boolean('fragments', gid_objs, gid_tools,
                               remove_tool=remove_tool,
                               remove_obj=remove_obj,
                               keep_highest=keep_highest)    

    def boolean_union2d(self, input_entity, tool_entity,
                        removeObject=False, removeTool=False, delete=False):

        def get_dimtag(entity):
            dimtags = []
            for en in entity:
                dimtags.append(id2dimtag(en))
            return dimtags

        all_entity = input_entity + tool_entity
        out_entity = self._boolean_xxx('fuse', input_entity, tool_entity,
                                       removeObject=removeObject, removeTool=removeTool,
                                       delete=delete)

        self.factory.synchronize()

        out_dimtag = get_dimtag(out_entity)
        xmin, ymin, zmin, xmax, ymax, zmax = find_combined_bbox(
            self.model, out_dimtag)

        dprint1("bounding box", xmin, ymin, zmin, xmax, ymax, zmax)

        dx = xmax - xmin
        dy = ymax - ymin
        bbx = self.factory.addRectangle(xmin - dx / 10., ymin - dy / 10., (zmin + zmax) / 2.,
                                        dx * 1.2, dy * 1.2)
        out_dimtag2, dimtagMap = self.factory.cut(((2, bbx),), out_dimtag)
        # print(out_dimtag2)
        bbx = self.factory.addRectangle(xmin - dx / 10., ymin - dy / 10., (zmin + zmax) / 2.,
                                        dx * 1.2, dy * 1.2)
        out_dimtag3, dimtagMap = self.factory.cut(((2, bbx),), out_dimtag2)
        self.factory.synchronize()
        return dimtag2id(out_dimtag3)

    def apply_fragments(self):
        if len(self.solids) > 1:
            keys = self.solids.keys()
        elif len(self.solids) == 1:
            keys = []
        elif len(self.faces) > 1:
            keys = self.faces.keys()
        elif len(self.faces) == 1:
            keys = []
        elif len(self.edges) > 1:
            keys = self.edges.keys()
        else:
            keys = []

        if len(keys) > 1:
            gid_objs = keys[:1]
            gid_tools = keys[1:]
            
            self.fragments(gid_objs, gid_tools,
                         remove_obj=True, remove_tool=True)

    def remove(self, gid, recursive=True):
        topolist = self.get_topo_list_for_gid(gid)
        #topo_list_child = self.get_topo_list_for_gid(gid, child=1)

        shape = topolist[gid]
        children = list(topolist.get_children(gid))

        self.builder.Remove(self.shape, shape)

        if not recursive:
            child_mapper = topolist.get_chilld_mapper(self.shape)
            for child in children:
                flag = child_mapper.Contains(shape)
                if not flag:  # need to put it back
                    self.builder.Add(self.shape, child)

    def inverse_remove(self, gids):
        comp = TopoDS_Compound()
        b = self.builder
        b.MakeCompound(comp)

        for gid in gids:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            self.builder.Add(comp, shape)

        self.shape = comp

        return gids

    def copy(self, gid):
        topolist = self.get_topo_list_for_gid(gid)
        shape = topolist[gid]

        copier = BRepBuilderAPI_Copy()
        copier.Perform(shape)
        if not copier.IsDone():
            assert False, "Can not copy shape"
        shape = copier.Shape()

        self.builder.Add(self.shape, shape)

        gid = topolist.add(shape)
        return gid

    def _perform_transform(self, gid, transformer, copy):
        topolist = self.get_topo_list_for_gid(gid)
        shape = topolist[gid]

        transformer.Perform(shape, copy)

        if not transformer.IsDone():
            assert False, "can not translate"

        new_shape = transformer.ModifiedShape(shape)
        isNew = not new_shape.IsSame(shape)

        if isNew:
            self.builder.Add(self.shape, new_shape)
            new_gid = topolist.add(new_shape)
        else:
            new_gid = None

        if not copy and isNew:
            #  I am not sure why OCC always copy object
            self.remove(gid)

        return new_gid

    def translate(self, gid, delta, copy=False):
        trans = gp_Trsf()
        trans.SetTranslation(gp_Vec(delta[0], delta[1], delta[2]))
        transformer = BRepBuilderAPI_Transform(trans)

        return self._perform_transform(gid, transformer, copy)

    def rotate(self, gid, point_on_axis, axis_dir, angle, copy=False):
        trans = gp_Trsf()

        x, y, z = point_on_axis
        ax, ay, az = axis_dir
        axis_revolution = gp_Ax1(gp_Pnt(x, y, z), gp_Dir(ax, ay, az))
        trans.SetRotation(axis_revolution, angle)

        transformer = BRepBuilderAPI_Transform(trans)

        return self._perform_transform(gid, transformer, copy)

    def dilate(self, gid, xyz, abc, copy=False):
        x, y, z = xyz
        a, b, c = abc
        gt = gp_GTrsf()
        gt.SetVectorialPart(gp_Mat(a, 0, 0, 0, b, 0, 0, 0, c))
        gt.SetTranslationPart(gp_XYZ(x * (1 - a), y * (1 - b), z * (1 - c)))

        transformer = BRepBuilderAPI_GTransform(gt)
        return self._perform_transform(gid, transformer, copy)

    def symmetrize(self, gid, abcd, copy=False):
        a, b, c, d = abcd
        gt = gp_GTrsf()
        p = max((a * a + b * b + c * c), 1e-12)
        f = -2.0 / p
        vec = (a * d * f, b * d * f, c * d * f)
        mat = (1 + a * a * f,
               a * b * f,
               a * c * f,
               a * b * f,
               1. + b * b * f,
               b * c * f,
               a * c * f,
               b * c * f,
               1. + c * c * f)
        gt.SetVectorialPart(gp_Mat(*mat))
        gt.SetTranslationPart(gp_XYZ(*vec))

        transformer = BRepBuilderAPI_GTransform(gt)
        return self._perform_transform(gid, transformer, copy)

    def extrude(self, gids, translation=None, rotation=None, wire=None):

        '''
        comp = TopoDS_Compound()
        self.builder.MakeCompound(comp);

        for gid in gids:
           topolist = self.get_topo_list_for_gid(gid)
           shape = topolist[gid]
           self.builder.Add(comp, shape)
        '''
        ret = []
        for gid in gids:
            topolist = self.get_topo_list_for_gid(gid)
            shape = topolist[gid]
            delete_input = topolist.is_toplevel(gid, self.shape)

            if translation is not None:
                dx, dy, dz = translation
                p = BRepPrimAPI_MakePrism(shape, gp_Vec(dx, dy, dz), False)

            elif rotation is not None:
                x, y, z = rotation[0]
                ax, ay, az = rotation[1]
                angle = rotation[2]
                pnt = gp_Pnt(x, y, z)
                dr = gp_Dir(ax, ay, az)
                ax = gp_Ax1(pnt, dr)
                p = BRepPrimAPI_MakeRevol(shape, ax, angle, False)

            elif wire is not None:
                from OCC.Core.GeomFill import GeomFill_IsDiscreteTrihedron
                p = BRepOffsetAPI_MakePipe(wire, shape, GeomFill_IsDiscreteTrihedron)

            else:
                assert False, "unknonw option"

            p.Build()
            if not p.IsDone():
                assert False, "can not extrude : " + gid

            if delete_input:
                self.builder.Remove(self.shape, shape)

            last = p.LastShape()

            if translation is not None:
                result = p.Prism().Shape()
            else:
                result = p.Shape()

            gid_last = self.add_to_topo_list(last)
            gid_extruded = self.add_to_topo_list(result)

            ret.append((gid_last, gid_extruded, result))

        return ret

    '''
    def update_topo_list_from_history(self, operator, list_of_shapes):
        iterator = TopTools_ListIteratorOfListOfShape(list_of_shapes)
        while iterator.More():
            shape = iterator.Value()            
            iterator.Next()
            # Do I need to do something with modified?
            # operator.Modified(shape)

            shape_gone = operator.IsDeleted(shape)
            if shape_gone:
                print("shape gone", shape_gone)

            shapes_new = operator.Generated(shape)
            iterator2 = TopTools_ListIteratorOfListOfShape(shapes_new)            
            while iterator2.More():
                shape_new = iterator2.Value()
                iterator2.Next()
                print("shape new", shape_new)
    '''
    def defeature(self, gid, gids_face):

        aSolid = self.solids[gid]
        
        features = TopTools_ListOfShape()
        for tmp in gids_face:
            topolist = self.get_topo_list_for_gid(tmp)
            shape = topolist[tmp]
            features.Append(shape)

        aDF = BRepAlgoAPI_Defeaturing()
        aDF.SetShape(aSolid)
        aDF.AddFacesToRemove(features)
        aDF.SetRunParallel(self.occ_parallel)
        aDF.SetToFillHistory(False)
        aDF.Build()

        if not aDF.IsDone():
            assert False, "Cannot remove faces"

        result = aDF.Shape()
        self.remove(gid)

        result = self.select_highest_dim(result)
        new_objs = self.register_shaps_balk(result)

        return new_objs

    def add_sequence(self, gui_name, gui_param, geom_name):
        self.geom_sequence.append((gui_name, gui_param, geom_name))


    '''
    high level interface:
       methods below directroy corresponds to GUI interface
       these routine should
        1) call builder.Add to add a shape to Compound
        2) register the name of shape
    '''


    ## 0D vertices
    def Point_build_geom(self, objs, *args):
        xarr, yarr, zarr = args

        _newobjs = []
        try:
            pos = np.vstack((xarr, yarr, zarr)).transpose()
        except:
            assert False, "can not make proper input array"

        PTs = [self.add_point(p) for p in pos]

        for p in PTs:
            shape = self.vertices[p]
            self.builder.Add(self.shape, shape)
            newkey = objs.addobj(p, 'pt')
            _newobjs.append(newkey)

        return list(objs), _newobjs

    ## 1D edges
    def Line_build_geom(self, objs, *args):
        xarr, yarr, zarr, make_spline, periodic = args
        lcar = 0.0
        if len(xarr) < 2:
            return
        try:
            pos = np.vstack((xarr, yarr, zarr)).transpose()
        except BaseException:
            assert False, "can not make proper input array"

        dist = np.sqrt(np.sum((pos[:-1, :] - pos[1:, :])**2, 1))

        if min(dist) == 0.0:
            assert False, "minimum distance between point is 0.0"
        if max(dist) > min(dist) * 1e4:
            assert False, "some points are too close (d_max > d_min*1e4)"

        if not make_spline:
            pts = [self.add_point(p) for ii, p in enumerate(pos)]
            pts1 = pts[:-1]
            pts2 = pts[1:]

            newkeys = []
            for p1, p2 in zip(pts1, pts2):
                ln = self.add_line(p1, p2)
                shape = self.edges[ln]
                self.builder.Add(self.shape, shape)
                newkeys.append(objs.addobj(ln, 'ln'))
            if periodic:
                ln = self.add_line(pts[-1], pts[0])
                shape = self.edges[ln]
                self.builder.Add(self.shape, shape)
                newkeys.append(objs.addobj(ln, 'ln'))

            _newobjs = newkeys
            if not periodic:
                newobj1 = objs.addobj(pts[0], 'pt')
                newobj2 = objs.addobj(pts[-1], 'pt')
                _newobjs.append(newobj1)
                _newobjs.append(newobj2)
        else:
            spline = self.add_spline(pos, periodic=periodic)
            shape = self.edges[spline]
            self.builder.Add(self.shape, shape)

            newobj = objs.addobj(spline, 'sp')
            _newobjs = [newobj]

        return list(objs), _newobjs

    def Polygon_build_geom(self, objs, *args):
        assert False, "polygon is not available in OCC geometry wrapper"

    def Spline_build_geom(self, objs, *args):
        pts = args
        pts = [x.strip() for x in pts[0].split(',')]
        gids = self.get_target1(objs, pts, 'p')

        if len(gids) < 3:
            assert False, "Spline requires more than 2 guide points"
        if len(gids) == 2 and gids[0] == gids[-1]:
            assert False, "Spline loop requires more than 3 guide points"
            
        if len(gids) > 3 and gids[0] == gids[-1]:
            periodic = True
            gids = gids[:-1]
        else:
            periodic = False        
            
        pos = np.vstack([self.get_point_coord(gid) for gid in gids])

        spline = self.add_spline(pos, periodic=periodic)
        shape = self.edges[spline]
        self.builder.Add(self.shape, shape)

        newobj = objs.addobj(spline, 'sp')
        return list(objs), [newobj]

    def CreateLine_build_geom(self, objs, *args):
        pts = args
        pts = [x.strip() for x in pts[0].split(',')]
        gids = self.get_target1(objs, pts, 'p')

        newkeys = []
        
        for i in range(len(gids)-1):
            p0 = gids[i]
            p1 = gids[i+1]            
            ln = self.add_line(p0, p1)
            shape = self.edges[ln]
            self.builder.Add(self.shape, shape)
            newkeys.append(objs.addobj(ln, 'ln'))

        return list(objs), newkeys

    def LineLoop_build_geom(self, objs, *args):
        assert False, "We don't support this"        

    ## 2D faces
    def Rect_build_geom(self, objs, *args):
        c1, e1, e2 = args
        lcar = 0.0

        c1 = np.array(c1)
        e1 = np.array(e1)
        e2 = np.array(e2)
        p1 = self.add_point(c1)
        p2 = self.add_point(c1 + e1)
        p3 = self.add_point(c1 + e1 + e2)
        p4 = self.add_point(c1 + e2)
        l1 = self.add_line(p1, p2)
        l2 = self.add_line(p2, p3)
        l3 = self.add_line(p3, p4)
        l4 = self.add_line(p4, p1)
        ll1 = self.add_line_loop([l1, l2, l3, l4])
        rec1 = self.add_plane_surface(ll1)

        shape = self.faces[rec1]
        self.builder.Add(self.shape, shape)

        newkey = objs.addobj(rec1, 'rec')
        return list(objs), [newkey]

    def Circle_build_geom(self, objs, *args):
        center, ax1, ax2, radius = args

        a1 = np.array(ax1)
        a2 = np.array(ax2)
        a2 = np.cross(np.cross(a1, a2), a1)
        a1 = a1 / np.sqrt(np.sum(a1**2)) * radius
        a2 = a2 / np.sqrt(np.sum(a2**2)) * radius

        c = np.array(center)
        p1 = self.add_point(c + a1)
        p2 = self.add_point(c + a2)
        p3 = self.add_point(c - a1)
        p4 = self.add_point(c - a2)
        ca1 = self.add_circle_arc(p1, p2, p3)
        ca2 = self.add_circle_arc(p3, p4, p1)
        ll1 = self.add_line_loop([ca1, ca2])

        ps1 = self.add_plane_surface(ll1)

        shape = self.faces[ps1]
        self.builder.Add(self.shape, shape)

        self.synchronize_topo_list(action='both')
        newkey = objs.addobj(ps1, 'ps')
        return list(objs), [newkey]

    def CreateSurface_build_geom(self, objs, *args):
        pts, isFilling = args
        pts = [x.strip() for x in pts.split(',')]

        gids_edge = self.get_target1(objs, pts, 'l')

        if isFilling:
            gids_vertex = []
            face_id = self.add_surface_filling(gids_edge, gids_vertex)            
            shape = self.faces[face_id]
            self.builder.Add(self.shape, shape)
            newobj2 = objs.addobj(face_id, 'sf')
            newkeys = [newobj2]
            '''
            gids_new = self.add_plate_surface(gids_edge, gids_vertex)

            newkeys = []
            for gid in gids_new:
                newkeys.append(objs.addobj(gid, 'sf'))
            '''

        else:
            ill = self.add_line_loop(edges)
            ips = self.add_plane_surface(ill)

            shape = self.faces[ips]
            self.builder.Add(self.shape, shape)
            newobj2 = objs.addobj(ips, 'ps')
            newkeys = [newobj2]
            
        return list(objs), newkeys

    def SurfaceLoop_build_geom(self, objs, *args):
        assert False, "We don't support this"

    ## 3D solids
    def CreateVolume_build_geom(self, objs, *args):

        pts = args
        pts = [x.strip() for x in pts[0].split(',')]

        gids = self.get_target1(objs, pts, 'f')
        sl = self.add_surface_loop(gids)
        v1 = self.add_volume(sl)

        shape = self.solids[v1]
        self.builder.Add(self.shape, shape)

        newobj2 = objs.addobj(v1, 'vol')

        return list(objs), [newobj2]


    def Box_build_geom(self, objs, *args):
        c1, e1, e2, e3 = args
        lcar = 0.0
        c1 = np.array(c1)
        e1 = np.array(e1)
        e2 = np.array(e2)
        p1 = c1
        p2 = c1 + e1
        p3 = c1 + e2
        p4 = c1 + e3
        p5 = c1 + e1 + e2
        p6 = c1 + e2 + e3
        p7 = c1 + e3 + e1
        p8 = c1 + e3 + e2 + e1

        v1 = self.add_box((p1, p2, p3, p4, p5, p6, p7, p8,))
        shape = self.solids[v1]
        self.builder.Add(self.shape, shape)

        newkey = objs.addobj(v1, 'bx')
        return list(objs), [newkey]

    def Ball_build_geom(self, objs, *args):

        x0, l1, l2, l3, a1, a2, a3 = args
        radii = [l1, l2, l3]
        rr = min(radii)

        gids_new = self.add_sphere(x0, rr, a1/180*np.pi, a2/180*np.pi, a3/180*np.pi)
        newkeys = []

        ss = (l1/rr, l2/rr, l3/rr)
        if ss[0] != ss[1] or ss[1] != ss[2]:
            gids_new = [self.dilate(gids_new[0], x0, ss, copy=False)]

        for gid_new in gids_new:
             newkeys.append(objs.addobj(gid_new, 'bl'))

        self.synchronize_topo_list(action='both')

        return list(objs), newkeys

    def Cone_build_geom(self, objs, *args):
        x0, d0, r1, r2, angle = args

        gids_new = self.add_cone(x0[0], x0[1], x0[2], d0[0], d0[1], d0[2],
                           r1, r2, angle/180*np.pi)

        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'cn'))

        return list(objs), newkeys

    def Cylinder_build_geom(self, objs, *args):
        x0, d0, r1, angle = args

        gids_new = self.add_cylinder(x0, d0, r1, angle/180*np.pi)
        
        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'cyl'))

        print(self.solids.d, gids_new)

        return list(objs), newkeys

        '''
        a1 = a1 / np.sqrt(np.sum(a1**2)) * r1
        a2 = a2 / np.sqrt(np.sum(a2**2)) * r1

        c = np.array(x0)
        p1 = self.add_point(c + a1, lcar)
        p2 = self.add_point(c + a2, lcar)
        p3 = self.add_point(c - a1, lcar)
        p4 = self.add_point(c - a2, lcar)
        pc = self.add_point(c, lcar)
        ca1 = self.add_circle_arc(p1, pc, p2)
        ca2 = self.add_circle_arc(p2, pc, p3)
        ca3 = self.add_circle_arc(p3, pc, p4)
        ca4 = self.add_circle_arc(p4, pc, p1)
        ll1 = self.add_line_loop([ca1, ca2, ca3, ca4])
        ps1 = self.add_plane_surface(ll1)

        ret = self.extrude(ps1, translation=d0,)
        '''
    def Wedge_build_geom(self, objs, *args):
        x0, d0, ltx = args
        gids_new = self.add_wedge(x0, d0, ltx)

        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'wg'))

        return list(objs), newkeys
    
    def Torus_build_geom(self, objs, *args):
        x0, r1, r2, angle, keep_interior = args

        gids_new = self.add_torus(x0, r1, r2, angle*np.pi/180)
        
        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'trs'))

        return list(objs), newkeys

    ## prutrusions
    def Extrude_build_geom(self, objs, *args):
        targets, tax, lengths = args
        
        targets = [x.strip() for x in targets.split(',')]
        gids = self.get_target2(objs, targets)

        print("tax", tax)

        trans = []
        if tax[0] == 'normal' or tax[0] == 'normalp':
            assert isinstance(gid[0], SurfaceID), "target must be surface"
            n1, p0 = self.get_planesurface_normal(gid[0])
             
            for length in lengths:
                if tax[1]:
                    tt = -n1 * length
                else:
                    tt = n1 * length
                trans.append(tt)
                
        elif tax[0] == 'normalp':
            assert len(lengths) == 1, "length should have one element"
            assert isinstance(gid[0], SurfaceID), "target must be surface"
            n1, p0 = self.get_planesurface_normal(gid[0])

            dests = [x.strip() for x in tax[1].split(',')]
            gid_dests = self.get_target2(objs, dests)

            for gid_dest in gid_dests:
                p1 = self.get_point_coord(gid_dest)
                if tax[2]:
                    tt = -n1 * np.sum((p1 - p0) * n1) * length
                else:
                    tt = n1 * np.sum((p1 - p0) * n1) * length
                trans.append(tt)

        elif tax[0] == 'fromto_points':
            dests1 = [x.strip() for x in tax[1].split(',')]
            dests2 = [x.strip() for x in tax[2].split(',')]

            gid_dests1 = self.get_target2(objs, dests1)
            gid_dests2 = self.get_target2(objs, dests2)

            assert len(gid_dests1) > 1, "Incorrect destination setting"
            assert len(gid_dests2) > 1, "Incorrect destination setting"

            p1 = self.get_point_coord(gid_dests1[0])
            p2 = self.get_point_coord(gid_dests2[0])

            n1 = p2 - p1
            if not tax[3]:
                n1 /= np.sqrt(np.sum(n1**2))
            if tax[4]:
                n1 *= -1
            for length in lengths:
                trans.append(length*n1)

        else:
            tax = np.array(tax).flatten()
            tax = tax / np.sqrt(np.sum(np.array(tax)**2))
            for length in lengths:
                trans.append(length*tax)

        newkeys = []

        for tt in trans:
            new_shapes = self.extrude(gids, translation=tt)

            gids = []
            for t, ret in zip(targets, new_shapes):
                gid_last, gid_extruded, shape = ret
                newkeys.append(objs.addobj(gid_last, t))
                newkeys.append(objs.addobj(gid_extruded, 'ex'))
                self.builder.Add(self.shape, shape)

                gids.append(gid_last)

        self.synchronize_topo_list(action='add')
        return list(objs), newkeys

    def Revolve_build_geom(self, objs, *args):

        targets, pax, rax, angles = args

        targets = [x.strip() for x in targets.split(',')]
        gids = self.get_target2(objs, targets)

        newkeys = []

        for angle in angles:
            rot = (pax, rax, angle * np.pi/180)
            new_shapes = self.extrude(gids, rotation=rot)

            gids = []
            for t, ret in zip(targets, new_shapes):
                gid_last, gid_extruded, shape = ret
                newkeys.append(objs.addobj(gid_last, t))
                newkeys.append(objs.addobj(gid_extruded, 'ex'))
                self.builder.Add(self.shape, shape)

                gids.append(gid_last)

        self.synchronize_topo_list(action='add')
        return list(objs), newkeys

    def Sweep_build_geom(self, objs, *args):
        print("objs", objs)
        targets, lines = args
        targets = [x.strip() for x in targets.split(',')]

        gids = self.get_target2(objs, targets)

        lines = [x.strip() for x in lines.split(',')]
        gid_lines = self.get_target1(objs, lines, 'l')

        wire_id = self.add_line_loop(gid_lines)
        wire = self.wires[wire_id]

        new_shapes = self.extrude(gids, wire=wire)
        newkeys = []

        for t, ret in zip(targets, new_shapes):
            gid_last, gid_extruded, shape = ret
            newkeys.append(objs.addobj(gid_last, t))
            newkeys.append(objs.addobj(gid_extruded, 'swp'))
            self.builder.Add(self.shape, shape)

        self.synchronize_topo_list(action='add')

        return list(objs), newkeys

    ## translation
    def Move_build_geom(self, objs, *args):
        targets, dx, dy, dz, keep = args
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        for gid in gids:
            new_gid = self.translate(gid, (dx, dy, dz), copy=keep)
            if new_gid is not None:
                newkeys.append(objs.addobj(new_gid, 'mv'))

        self.synchronize_topo_list(action='both')
        print('faces', self.faces.d)

        return list(objs), newkeys

    def Rotate_build_geom(self, objs, *args):
        targets, point_on_axis, axis_dir, angle, keep = args

        newkeys = []        
        targets = [x.strip() for x in targets.split(',')]

        gids = self.get_target2(objs, targets)

        for gid in gids:
            new_gid = self.rotate(gid, point_on_axis, axis_dir,
                                  np.pi * angle / 180., copy=keep)
            if new_gid is not None:
                newkeys.append(objs.addobj(new_gid, 'mv'))

        self.synchronize_topo_list(action='both')
        print('faces', self.faces.d)

        return list(objs), newkeys

    def Scale_build_geom(self, objs, *args):
        targets, cc, ss, keep = args
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        for gid in gids:
            new_gid = self.dilate(gid, cc, ss, copy=keep)
            if new_gid is not None:
                newkeys.append(objs.addobj(new_gid, 'sc'))

        self.synchronize_topo_list(action='both')

        return list(objs), newkeys

    def Flip_build_geom(self, objs, *args):
        targets, a, b, c, d, keep = args
        abcd = (a, b, c, d)
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        for gid in gids:
            new_gid = self.symmetrize(gid, abcd, copy=keep)
            if new_gid is not None:
                newkeys.append(objs.addobj(new_gid, 'flp'))

        self.synchronize_topo_list(action='both')

        return list(objs), newkeys

    def Array_build_geom(self, objs, *args):
        targets, count, displacement = args
        dx, dy, dz = displacement
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        i = 1
        while i < count:
            for gid in gids:
                delta = (dx * i, dy * i, dz * i)
                new_gid = self.translate(gid, delta, True)
                if new_gid is not None:
                    newkeys.append(objs.addobj(new_gid, 'cp'))
            i = i + 1

        self.synchronize_topo_list(action='add')

        return list(objs), newkeys

    def ArrayRot_build_geom(self, objs, *args):
        targets, count, point_on_axis, axis_dir, angle = args

        newkeys = []

        targets = [x.strip() for x in targets.split(',')]
        gids = self.get_target2(objs, targets)

        i = 1
        while i < count:
            for gid in gids:
                angle1 = angle * i
                new_gid = self.rotate(gid, point_on_axis, axis_dir,
                                      np.pi * angle1 / 180., True)

                if new_gid is not None:
                    newkeys.append(objs.addobj(new_gid, 'cp'))
            i = i + 1

        self.synchronize_topo_list(action='add')

        return list(objs), newkeys

    ## fillet/chamfer                           
    def Fillet_build_geom(self, objs, *args):

        volumes, curves, radii = args
        volumes = [x.strip() for x in volumes.split(',')]
        curves = [x.strip() for x in curves.split(',')]

        gid_vols = self.get_target1(objs, volumes, 'v')
        gid_curves = self.get_target1(objs, curves, 'l')

        gids_new = self.fillet(gid_vols, gid_curves, radii)

        newkeys = []        
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'vol'))

        self.synchronize_topo_list()
        return list(objs), newkeys

    def Chamfer_build_geom(self, objs, *args):
        volumes, curves, distances, surfaces = args

        volumes = [x.strip() for x in volumes.split(',')]
        curves = [x.strip() for x in curves.split(',')]
        surfaces = [x.strip() for x in surfaces.split(',')]

        gid_vols = self.get_target1(objs, volumes, 'v')
        gid_curves = self.get_target1(objs, curves, 'l')
        gid_faces = self.get_target1(objs, surfaces, 'f')

        gids_new = self.chamfer(gid_vols, gid_curves, gid_faces, distances)

        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'vol'))

        self.synchronize_topo_list()
        return list(objs), newkeys

    ## copy/remove                           
    def Copy_build_geom(self, objs, *args):
        targets = args[0]
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)
        
        for gid in gids:
            copied_gid = self.copy(gid)
            newkeys.append(objs.addobj(copied_gid, 'cp'))

        self.synchronize_topo_list(action='add')
        return list(objs), newkeys

    def Remove_build_geom(self, objs, *args):
        targets, recursive = args
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)
        if len(gids) == 0:
            assert False, "empty imput objects: "+','.join(targets) 

        for gid in gids:
            self.remove(gid, recursive=recursive)
        self.synchronize_topo_list()

        for t in targets:
            if t in objs:
                del objs[t]

        return list(objs), newkeys

    def Remove2_build_geom(self, objs, *args):
        targets = args[0]
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        ret = self.inverse_remove(gids)
        self.synchronize_topo_list()

        for t in list(objs):
            del objs[t]
        for rr in ret:
            newkeys.append(objs.addobj(rr, 'kpt'))

        return list(objs), newkeys
                               
    def RemoveFaces_build_geom(self, objs, *args):
        targets, faces = args
        
        targets = [x.strip() for x in targets.split(',')]
        faces = [x.strip() for x in faces.split(',')]

        if len(targets) != 1:
            assert False, "Chose one volume"

        gid = self.get_target1(objs, targets, 'v')[0]
        gids_face = self.get_target1(objs, faces, 'f')

        gids_new = self.defeature(gid, gids_face)

        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'dftr'))

        self.synchronize_topo_list(verbose=True)
        return list(objs), newkeys

    def Union_build_geom(self, objs, *args):
        tp, tm, delete_input, delete_tool, keep_highest = args
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]

        gid_objs = self.get_target2(objs, tp)
        gid_tools = self.get_target2(objs, tm)

        gids_new = self.union(gid_objs, gid_tools,
                              remove_obj=delete_input,
                              remove_tool=delete_tool,
                              keep_highest=keep_highest)

        newkeys = []
        for gid in gids_new:
            newkeys.append(objs.addobj(gid, 'uni'))

        self.synchronize_topo_list(verbose=True)
        return list(objs), newkeys


    def Difference_build_geom(self, objs, *args):
        tp, tm, delete_input, delete_tool, keep_highest = args
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]

        gid_objs = self.get_target2(objs, tp)
        gid_tools = self.get_target2(objs, tm)

        gids_new = self.difference(gid_objs, gid_tools,
                                   remove_obj=delete_input,
                                   remove_tool=delete_tool,
                                   keep_highest=keep_highest)

        newkeys = []
        for gid in gids_new:
            #topolist = self.get_topo_list_for_gid(gid)
            #shape = topolist[gid]
            newkeys.append(objs.addobj(gid, 'diff'))

        self.synchronize_topo_list(verbose=True)
        return list(objs), newkeys

    def Intersection_build_geom(self, objs, *args):
        tp, tm, delete_input, delete_tool, keep_highest = args
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]

        gid_objs = self.get_target2(objs, tp)
        gid_tools = self.get_target2(objs, tm)

        gids_new = self.intersection(gid_objs, gid_tools,
                                     remove_obj=delete_input,
                                     remove_tool=delete_tool,
                                     keep_highest=keep_highest)

        newkeys = []
        for gid in gids_new:
            #topolist = self.get_topo_list_for_gid(gid)
            #shape = topolist[gid]
            newkeys.append(objs.addobj(gid, 'diff'))

        self.synchronize_topo_list(verbose=True)
        return list(objs), newkeys

    def Fragments_build_geom(self, objs, *args):
        tp, tm, delete_input, delete_tool, keep_highest = args
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]

        gid_objs = self.get_target2(objs, tp)
        gid_tools = self.get_target2(objs, tm)

        gids_new = self.fragments(gid_objs, gid_tools,
                                  remove_obj=delete_input,
                                  remove_tool=delete_tool,
                                  keep_highest=keep_highest)

        newkeys = []
        for gid in gids_new:
            #topolist = self.get_topo_list_for_gid(gid)
            #shape = topolist[gid]
            newkeys.append(objs.addobj(gid, 'diff'))

        self.synchronize_topo_list(verbose=True)
        return list(objs), newkeys


        '''    
        input_entity = get_target2(objs, tp)
        tool_entity = get_target2(objs, tm)
        ret = self.boolean_fragments(
            input_entity,
            tool_entity,
            removeObject=delete_input,
            removeTool=delete_tool)

        newkeys = []
        newkeys = []
        for rr in ret:
            if rr.dim == input_entity[0].dim:
                newkeys.append(objs.addobj(rr, 'frag'))
            else:
                if keep_highest:
                    self.remove([rr], recursive=True)
                else:
                    newkeys.append(objs.addobj(rr, get_geom_key(rr)))

        if delete_input:
            for x in tp:
                if x in objs:
                    del objs[x]
        if delete_tool:
            for x in tm:
                if x in objs:
                    del objs[x]

        return list(objs), newkeys
        '''
    '''
    2D elements
    '''

    def Point2D_build_geom(self, objs, *args):
        xarr, yarr = args
        xarr = np.atleast_1d(xarr)
        yarr = np.atleast_1d(yarr)
        zarr = xarr * 0.0
        try:
            pos = np.vstack((xarr, yarr, zarr)).transpose()
        except:
            assert False, "can not make proper input array"

        PTs = [self.add_point(p) for p in pos]

        for p in PTs:
            shape = self.vertices[p]
            self.builder.Add(self.shape, shape)
            newkey = objs.addobj(p, 'pt')
            _newobjs.append(newkey)

        return list(objs), _newobjs

    # Define 2D version the same as 3D
    Line2D_build_geom = Line_build_geom

    def Circle2D_build_geom(self, objs, *args):
        center, ax1, ax2, radius = args

        a1 = np.array(ax1 + [0])
        a2 = np.array(ax2 + [0])
        a2 = np.cross(np.cross(a1, a2), a1)
        a1 = a1 / np.sqrt(np.sum(a1**2)) * radius
        a2 = a2 / np.sqrt(np.sum(a2**2)) * radius

        c = np.array(center + [0])
        p1 = self.add_point(c + a1)
        p2 = self.add_point(c + a2)
        p3 = self.add_point(c - a1)
        p4 = self.add_point(c - a2)
        ca1 = self.add_circle_arc(p1, p2, p3)
        ca2 = self.add_circle_arc(p3, p4, p1)        
        ll1 = self.add_line_loop([ca1, ca2])

        ps1 = self.add_plane_surface(ll1)

        shape = self.faces[ps1]
        self.builder.Add(self.shape, shape)

        self.synchronize_topo_list(action='both')
        newkey = objs.addobj(ps1, 'ps')

        return list(objs), [newkey]

    def Arc2D_build_geom(self, objs, *args):
        center, ax1, ax2, radius, an1, an2, do_fill = args
        lcar = 0.0
        a1 = np.array(ax1 + [0])
        a2 = np.array(ax2 + [0])
        a2 = np.cross(np.cross(a1, a2), a1)
        a1 = a1 / np.sqrt(np.sum(a1**2)) * radius
        a2 = a2 / np.sqrt(np.sum(a2**2)) * radius
        if an1 > an2:
            tmp = an2
            an2 = an1
            an1 = tmp

        if an2 - an1 >= 360:
            assert False, "angle must be less than 360"

        an3 = (an1 + an2) / 2.0
        pt1 = a1 * np.cos(an1 * np.pi / 180.) + a2 * np.sin(an1 * np.pi / 180.)
        pt2 = a1 * np.cos(an2 * np.pi / 180.) + a2 * np.sin(an2 * np.pi / 180.)
        pt3 = a1 * np.cos(an3 * np.pi / 180.) + a2 * np.sin(an3 * np.pi / 180.)

        c = np.array(center + [0])
        p1 = self.add_point(c + pt1, lcar)
        p2 = self.add_point(c + pt2, lcar)
        p3 = self.add_point(c + pt3, lcar)
        pc = self.add_point(c, lcar)
        ca1 = self.add_circle_arc(p1, pc, p3)
        ca2 = self.add_circle_arc(p3, pc, p2)

        if not do_fill:
            newkey1 = objs.addobj(ca1, 'ln')
            newkey2 = objs.addobj(ca2, 'ln')
            newkeys = [newkey1, newkey2]

        else:
            l1 = self.add_line(pc, p1)
            l2 = self.add_line(p2, pc)
            ll1 = self.add_line_loop([l1, ca1, ca2, l2])
            ps1 = self.add_plane_surface(ll1)
            newkeys = [objs.addobj(ps1, 'ps')]

        return list(objs), newkeys

    def Rect2D_build_geom(self, objs, *args):
        c1, e1, e2 = args
        lcar = 0.0
        c1 = np.array(c1 + [0])
        e1 = np.array(e1 + [0])
        e2 = np.array(e2 + [0])
        p1 = self.add_point(c1)
        p2 = self.add_point(c1 + e1)
        p3 = self.add_point(c1 + e1 + e2)
        p4 = self.add_point(c1 + e2)
        l1 = self.add_line(p1, p2)
        l2 = self.add_line(p2, p3)
        l3 = self.add_line(p3, p4)
        l4 = self.add_line(p4, p1)
        ll1 = self.add_line_loop([l1, l2, l3, l4])
        rec1 = self.add_plane_surface(ll1)

        shape = self.faces[rec1]
        self.builder.Add(self.shape, shape)

        newkey = objs.addobj(rec1, 'rec')
        return list(objs), [newkey]

    def Polygon2D_build_geom(self, objs, *args):
        del objs
        del args
        assert False, "We dont support this"

    def Move2D_build_geom(self, objs, *args):
        targets, dx, dy, keep = args
        dz = 0.0

        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        gids = self.get_target2(objs, targets)

        for gid in gids:
            new_gid = self.translate(gid, (dx, dy, dz), keep)
            if new_gid is not None:
                newkeys.append(objs.addobj(new_gid, 'mv'))

        self.synchronize_topo_list(action='both')

        return list(objs), newkeys

    def Rotate2D_build_geom(self, objs, *args):
        targets, cc, angle, keep = args
        cx, cy = cc
        cz = 0.0
        ax, ay, az = 0.0, 0.0, 1.0
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = get_target2(objs, targets)

        if keep:
            tt = self.copy(tt)
        self.rotate(tt, cx, cy, cz, ax, ay, az, np.pi * angle / 180.)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'rot'))

        return list(objs), newkeys

    def Flip2D_build_geom(self, objs, *args):
        targets, a, b, d, keep = args
        c = 0.0
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = get_target2(objs, targets)
        if keep:
            tt = self.copy(tt)
        self.symmetrize(tt, a, b, c, d)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'flp'))

        return list(objs), newkeys

    def Scale2D_build_geom(self, objs, *args):
        targets, cc, ss, keep = args
        cx, cy = cc
        cz = 0.0
        sx, sy = ss
        sz = 1.0
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = get_target2(objs, targets)
        if keep:
            tt = self.copy(tt)
        self.dilate(tt, cx, cy, cz, sx, sy, sz)
        if keep:
            for t in tt:
                newkeys.append(objs.addobj(t, 'sc'))

        return list(objs), newkeys

    def Array2D_build_geom(self, objs, *args):
        targets, count, displacement = args
        dx, dy = displacement
        dz = 0.0
        targets = [x.strip() for x in targets.split(',')]

        newkeys = []
        tt = get_target2(objs, targets)
        for i in range(count):
            tt = self.copy(tt)
            self.translate(tt, dx, dy, dz)
            for t in tt:
                newkeys.append(objs.addobj(t, 'cp'))

        return list(objs), newkeys

    def Union2D_build_geom(self, objs, *args):
        tp, tm, delete_input, delete_tool, keep_highest = args
        tp = [x.strip() for x in tp.split(',')]
        tm = [x.strip() for x in tm.split(',')]

        input_entity = get_target2(objs, tp)
        tool_entity = get_target2(objs, tm)

        ret = self.boolean_union2d(
            input_entity,
            tool_entity,
            removeObject=delete_input,
            removeTool=delete_tool)

        newkeys = []
        for rr in ret:
            if rr.dim == input_entity[0].dim:
                newkeys.append(objs.addobj(rr, 'uni'))
            else:
                if keep_highest:
                    self.remove([rr], recursive=True)
                else:
                    newkeys.append(objs.addobj(rr, get_geom_key(rr)))

        if delete_input:
            for x in tp:
                if x in objs:
                    del objs[x]

        if delete_tool:
            for x in tm:
                if x in objs:
                    del objs[x]

        return list(objs), newkeys

    def SplitByPlane_build_geom(self, objs, *args):
        print(args)

        def project_ptx_2_plain(normal, cptx, p):
            dp = p - cptx
            dp = dp - np.sum(dp * normal) * normal
            return dp + cptx

        def containing_bbox(normal, cptx, xmin, ymin, zmin, xmax, ymax, zmax):
            corners = (np.array([xmin, ymin, zmin]),
                       np.array([xmin, ymin, zmax]),
                       np.array([xmin, ymax, zmin]),
                       np.array([xmax, ymin, zmin]),
                       np.array([xmax, ymax, zmin]),
                       np.array([xmin, ymax, zmax]),
                       np.array([xmax, ymin, zmax]),
                       np.array([xmax, ymax, zmax]),)
            # projected point
            p = [project_ptx_2_plain(normal, cptx, pp) for pp in corners]

            # distance from plain
            d = [np.sum((pp - cptx) * normal) for pp in corners]
            dist1 = np.max(d)

            # distance on the plain
            d = [np.sqrt(np.sum((pp - cptx)**2)) for pp in p]
            idx = np.argmax(d)
            dist2 = np.max(d)

            size = np.max((dist1, dist2)) * 1.2

            n1 = (p[idx] - cptx)
            n2 = np.cross(normal, n1)

            c1 = cptx - n1 * size - n2 * size
            e1 = 2 * n1 * size
            e2 = 2 * n2 * size
            e3 = normal * size
            box = (c1, c1 + e1, c1 + e2, c1 + e3, c1 + e1 + e2, c1 + e2 + e3, c1 + e3 + e1,
                   c1 + e3 + e2 + e1)

            return box

        targets = [x.strip() for x in args[0].split(',')]
        tt = get_target2(objs, targets)

        dimtags = [id2dimtag(i) for i in tt]
        xmin, ymin, zmin, xmax, ymax, zmax = find_combined_bbox(
            self.model, dimtags)

        if args[1][0] == '3_points':
            # args[1] = ['3_points', '1', '7', '8']

            ptx1ID = get_target1(objs, [args[1][1], ], 'p')[0]
            ptx2ID = get_target1(objs, [args[1][2], ], 'p')[0]
            ptx3ID = get_target1(objs, [args[1][3], ], 'p')[0]
            ptx1 = np.array(gmsh.model.getValue(0, int(ptx1ID), []))
            ptx2 = np.array(gmsh.model.getValue(0, int(ptx2ID), []))
            ptx3 = np.array(gmsh.model.getValue(0, int(ptx3ID), []))

            n = np.cross(ptx1 - ptx2, ptx1 - ptx3)
            if np.sum(n**2) == 0:
                assert False, "three points does not span a surface."
            normal = n / np.sqrt(np.sum(n**2))
            cptx = (ptx1 + ptx2 + ptx3) / 3.0
        elif args[1][0] == 'by_abc':
            data = np.array(args[1][1]).flatten()
            xmin, ymin, zmin, xmax, ymax, zmax = find_combined_bbox(
                self.model, dimtags)
            normal = data[:3]
            xx = np.array(
                [(xmin + xmax) / 2, (ymin + ymax) / 2.0, (zmin + zmax) / 2.0])
            s = data[-1] - np.sum(normal * xx)
            cptx = xx + s * normal
        elif args[1][0] == 'face_parallel':
            faceID = get_target1(objs, [args[1][1], ], 'f')[0]
            ptxID = get_target1(objs, [args[1][2], ], 'p')[0]
            cptx = np.array(gmsh.model.getValue(0, int(ptxID), []))

            n1 = np.array(gmsh.model.getNormal(faceID, (0, 0)))
            n2 = np.array(gmsh.model.getNormal(faceID, (0, 1)))
            n3 = np.array(gmsh.model.getNormal(faceID, (1, 0)))
            n1 /= np.sqrt(np.sum(n1**2))
            n2 /= np.sqrt(np.sum(n2**2))
            n3 /= np.sqrt(np.sum(n3**2))

            if np.any(n1 != n2) or np.any(n1 != n3):
                assert False, "surface is not flat"
            normal = n1

        else:
            assert False, "unknown option:" + args

        points = containing_bbox(
            normal, cptx, xmin, ymin, zmin, xmax, ymax, zmax)
        v = self.add_box(points)

        ret1 = self.boolean_difference(tt, (v,),
                                       removeObject=False,
                                       removeTool=False)
        ret2 = self.boolean_intersection(tt, (v,),
                                         removeObject=True,
                                         removeTool=True)
        newkeys = []
        for rr in ret1 + ret2:
            if rr.dim == tt[0].dim:
                newkeys.append(objs.addobj(rr, 'splt'))
            else:
                if keep_highest:
                    self.remove([rr], recursive=True)
                else:
                    newkeys.append(objs.addobj(rr, get_geom_key(rr)))

        for x in targets:
            if x in objs:
                del objs[x]

        return list(objs), newkeys

        return list(objs), []

    def _WorkPlane_build_geom(self, objs, c1, a1, a2):
        x1 = np.array([1., 0., 0.])

        ax = np.cross(x1, a1)
        an = np.arctan2(np.sqrt(np.sum(ax**2)), np.dot(a1, x1))

        tt = [objs[t] for t in objs]

        #from petram.geom.gmsh_geom_wrapper import VertexID, LineID, SurfaceID

        tt = self.get_unique_entity(tt)

        #print("first rot ???", ax, an, np.sum(ax**2))
        if np.sum(ax**2) == 0.0:
            if an != 0.0:
                # if a1 is [0, 0, -1], rotate 180 deg
                ax = np.array([0, 1, 0])
                an = np.pi
            else:
                ax = x1
                an = 0.0
        if np.sum(ax**2) != 0.0 and an != 0.0:
            #print("first rot", ax, an)
            self.rotate(tt, 0, 0, 0, ax[0], ax[1], ax[2], an)

        from petram.geom.geom_utils import rotation_mat
        R = rotation_mat(ax, an)
        '''
        c = np.cos(an); s = np.sin(an)
        R = np.array(
            [[c + (1-c)*ax[0]**2, ax[0]*ax[1]*(1-c)-ax[2]*s, ax[0]*ax[2]*(1-c)+ax[1]*s],
             [ax[0]*ax[1]*(1-c)+ax[2]*s, c + (1-c)*ax[1]**2,  ax[1]*ax[2]*(1-c)-ax[0]*s],
             [ax[0]*ax[2]*(1-c)-ax[1]*s, ax[1]*ax[2]*(1-c)+ax[0]*s, c + (1-c)*ax[2]**2]]
            )
        '''
        y2 = np.dot(R, np.array([0, 1, 0]))
        ax = a1
        aaa = np.cross(a1, y2)
        an = np.arctan2(np.dot(a2, aaa), np.dot(a2, y2))

        # for t in tt:
        #     if isinstance(t, SurfaceID): continue
        #     print("working on t", t)
        #
        #     geom.rotate([t], 0, 0, 0, ax[0], ax[1], ax[2], an)
        #print("2nd rot ???", ax, an, np.sum(ax**2))
        if np.sum(ax**2) == 0.0 and an != 0.0:
            # rotate 180 deg around a1
            ax = a1
            an = np.pi
        if np.sum(ax**2) != 0.0 and an != 0.0:
            #print("2nd rot", ax, an)
            self.rotate(tt, 0, 0, 0, ax[0], ax[1], ax[2], an)

        if c1[0] != 0.0 or c1[1] != 0.0 or c1[2] != 0.0:
            self.translate(tt, c1[0], c1[1], c1[2])

        #self._newobjs = objs.keys()
        return list(objs), []

    def WorkPlane_build_geom(self, objs, *args):
        c1, a1, a2 = args
        c1 = np.array(c1)
        a1 = np.array(a1)
        a1 = a1 / np.sqrt(np.sum(a1**2))
        a2 = np.array(a2)
        a2 = a2 / np.sqrt(np.sum(a2**2))
        return self._WorkPlane_build_geom(objs, c1, a1, a2)

    def WorkPlaneByPoints_build_geom(self, objs, *args):
        c1, a1, a2, flip1, flip2 = args

        self.factory.synchronize()
        c1 = gmsh.model.getValue(0, int(c1), [])
        a1 = gmsh.model.getValue(0, int(a1), [])
        a2 = gmsh.model.getValue(0, int(a2), [])

        d1 = np.array(a1) - np.array(c1)
        d1 = d1 / np.sqrt(np.sum(d1**2))
        if flip1:
            d1 = -d1

        d2 = np.array(a2) - np.array(c1)
        d2 = d2 / np.sqrt(np.sum(d2**2))

        d3 = np.cross(d1, d2)
        d3 = d3 / np.sqrt(np.sum(d3**2))
        d2 = np.cross(d3, d1)
        d2 = d2 / np.sqrt(np.sum(d2**2))
        if flip2:
            d2 = -d2

        return self._WorkPlane_build_geom(objs, c1, d1, d2)

    def healShapes(self, dimtags, fix_tol, fixDegenerated=False,
                   fixSmallEdges=False,
                   fixSmallFaces=False,
                   sewFaces=False):

        self.factory.synchronize()
        top_level = self.get_toplevel_enteties()
        if dimtags is None:
            dimtags = top_level

        ret = []
        removed = []

        for dimtag in dimtags:
            if not dimtag in top_level:
                print(
                    "skipping " +
                    str(dimtag) +
                    " since it is not top level entitiy")
                continue
            outdimtags = self.factory.healShapes(dimTags=[dimtag],
                                                 tolerance=fix_tol,
                                                 fixDegenerated=fixDegenerated,
                                                 fixSmallEdges=fixSmallEdges,
                                                 fixSmallFaces=fixSmallFaces,
                                                 sewFaces=sewFaces)
            #print("heal outdimtags", outdimtags)
            self.factory.synchronize()
            self.factory.remove([dimtag], recursive=True)
            ret.append(outdimtags[0])
            removed.append(dimtag)

        self.factory.synchronize()
        return ret, removed

    def healCAD_build_geom(self, objs, *args):
        targets, use_fix_param, use_fix_tol = args

        self.factory.synchronize()

        targets = [x.strip()
                   for x in targets.split(',') if len(x.strip()) != 0]
        if len(targets) == 0:
            dimtags = None
        else:
            targetID = get_target2(objs, targets)
            dimtags = get_dimtag(targetID)

        ret, removed = self.healShapes(dimtags, use_fix_tol, fixDegenerated=use_fix_param[0],
                                       fixSmallEdges=use_fix_param[1],
                                       fixSmallFaces=use_fix_param[2],
                                       sewFaces=use_fix_param[3])

        for k in list(objs):
            if objs[k].to_dimtag() in removed:
                del objs[k]

        newkeys = []
        ret = dimtag2id(ret)
        for rr in ret:
            newkeys.append(objs.addobj(rr, 'hld'))

        return list(objs), newkeys

    def select_highest_dim(self, shape):
        comp = TopoDS_Compound()
        b = self.builder
        b.MakeCompound(comp)

        mmm = TopTools_IndexedMapOfShape()
        topexp_MapShapes(shape, TopAbs_SOLID, mmm)
        if mmm.Size() == 0:
            topexp_MapShapes(shape, TopAbs_FACE, mmm)
            if mmm.Size() == 0:
                topexp_MapShapes(shape, TopAbs_EDGE, mmm)
                if mmm.Size() == 0:
                    topexp_MapShapes(shape, TopAbs_VERTEX, mmm)
                    ex1 = TopExp_Explorer(shape, TopAbs_VERTEX)
                else:
                    ex1 = TopExp_Explorer(shape, TopAbs_EDGE)
            else:
                ex1 = TopExp_Explorer(shape, TopAbs_FACE)
        else:
            ex1 = TopExp_Explorer(shape, TopAbs_SOLID)

        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        return comp

    def register_shaps_balk(self, shape):
        solidMap, shellMap, faceMap, wireMap, edgeMap, vertMap = self.prep_maps(
            shape)

        '''
        solidMap = TopTools_IndexedMapOfShape()
        shellMap = TopTools_IndexedMapOfShape()
        faceMap = TopTools_IndexedMapOfShape()
        wireMap = TopTools_IndexedMapOfShape()
        edgeMap = TopTools_IndexedMapOfShape()
        vertMap = TopTools_IndexedMapOfShape()

        topexp_MapShapes(shape,TopAbs_SOLID, solidMap)
        topexp_MapShapes(shape,TopAbs_SHELL, shellMap)
        topexp_MapShapes(shape,TopAbs_FACE, faceMap)
        topexp_MapShapes(shape,TopAbs_WIRE, wireMap)
        topexp_MapShapes(shape,TopAbs_EDGE, edgeMap)
        topexp_MapShapes(shape,TopAbs_VERTEX,vertMap)
        '''

        usolids = topo_seen(mapping=solidMap)
        ushells = topo_seen(mapping=shellMap)
        ufaces = topo_seen(mapping=faceMap)
        uwires = topo_seen(mapping=wireMap)
        uedges = topo_seen(mapping=edgeMap)
        uvertices = topo_seen(mapping=vertMap)

        new_objs = []
        # registor solid
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

        register_topo(
            shape,
            ushells,
            TopAbs_SHELL,
            TopAbs_SOLID,
            topods_Shell,
            topods_Solid,
            self.shells)
        register_topo(
            shape,
            ufaces,
            TopAbs_FACE,
            TopAbs_SHELL,
            topods_Face,
            topods_Shell,
            self.faces,
            dim=2)
        register_topo(
            shape,
            uwires,
            TopAbs_WIRE,
            TopAbs_FACE,
            topods_Wire,
            topods_Face,
            self.wires)
        register_topo(
            shape,
            uedges,
            TopAbs_EDGE,
            TopAbs_WIRE,
            topods_Edge,
            topods_Wire,
            self.edges,
            dim=1)
        register_topo(shape, uvertices, TopAbs_VERTEX, TopAbs_EDGE, topods_Vertex, topods_Edge,
                      self.vertices, dim=0)

        b = self.builder
        comp = self.shape
        ex1 = TopExp_Explorer(shape, TopAbs_SOLID)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(shape, TopAbs_SHELL, TopAbs_SOLID)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(shape, TopAbs_FACE, TopAbs_SHELL)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(shape, TopAbs_WIRE, TopAbs_FACE)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(shape, TopAbs_EDGE, TopAbs_WIRE)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()
        ex1 = TopExp_Explorer(shape, TopAbs_VERTEX, TopAbs_EDGE)
        while ex1.More():
            b.Add(comp, ex1.Current())
            ex1.Next()

        return new_objs

    def importShape_common(self, shape, highestDimOnly,
                           use_fix_param, use_fix_tol, objs):
        '''
        (I need to address healing)
        if use_fix:
             PTs, void = self.healShapes(PTs, use_fix_tol, fixDegenerated = use_fix_param[0],
                                                fixSmallEdges = use_fix_param[1],
                                                fixSmallFaces = use_fix_param[2],
                                                sewFaces = use_fix_param[3])
        '''

        if highestDimOnly:
            shape = self.select_highest_dim(shape)
        new_objs = self.register_shaps_balk(shape)

        newkeys = []
        dim = max([p.idx for p in new_objs])
        for p in new_objs:
            if p.idx == dim:
                newkeys.append(objs.addobj(p, 'impt'))

        return list(objs), newkeys

    def BrepImport_build_geom(self, objs, *args):
        cad_file, use_fix, use_fix_param, use_fix_tol, highestDimOnly = args

        from OCC.Core.BRepTools import breptools_Read

        shape = TopoDS_Shape()
        success = breptools_Read(shape, cad_file, self.builder)

        if not success:
            assert False, "Failed to read brep"

        return self.importShape_common(
            shape, highestDimOnly, use_fix_param, use_fix_tol, objs)

    def CADImport_build_geom(self, objs, *args):
        from OCC.Core.STEPControl import STEPControl_Reader
        from OCC.Core.IGESControl import IGESControl_Reader
        from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
        from OCC.Core.Interface import Interface_Static_SetCVal

        unit = args[-1]
        cad_file, use_fix, use_fix_param, use_fix_tol, highestDimOnly = args[:-1]

        if (cad_file.lower().endswith(".iges") or
                cad_file.lower().endswith(".igs")):
            reader = IGESControl_Reader()
        elif (cad_file.lower().endswith(".step") or
              cad_file.lower().endswith(".stp")):
            reader = STEPControl_Reader()
        else:
            assert False, "unsupported format"

        check = Interface_Static_SetCVal("xstep.cascade.unit", unit)
        if not check:
            assert False, "can not set unit"

        status = reader.ReadFile(cad_file)

        if status == IFSelect_RetDone:  # check status
            failsonly = False
            reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
            reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
            reader.NbRootsForTransfer()
            reader.TransferRoot(1)
            shape = reader.Shape(1)
        else:
            assert False, "Error: can't read STEP file."

        return self.importShape_common(
            shape, highestDimOnly, use_fix_param, use_fix_tol, objs)

    def find_tinyloop(self, esize, thr=1e-5):
        # magic number to define too small
        el_max = max(list(esize.values()))
        edges = np.array([x for x in esize if esize[x] < el_max * thr])

        dimtags = self.model.getEntities(1)
        loops = []
        for dim, tag in dimtags:
            if len(self.model.getBoundary([(dim, tag)], oriented=False)) == 0:
                loops.append(tag)
        edges = [e for e in edges if e in loops]
        print("tiny loop edges ", edges)

        return edges

    def make_safe_file(self, filename, trash, ext):
        #map = self.getEntityNumberingInfo()
        # make filename safe
        filename = '_'.join(filename.split("/"))
        filename = '_'.join(filename.split(":"))
        filename = '_'.join(filename.split("\\"))

        if trash == '':  # when finalizing
            return os.path.join(os.getcwd(), filename + ext)
        else:
            return os.path.join(trash, filename + ext)

    def generate_preview_mesh(self, filename, trash, mesh_quality=1):
        if self.queue is not None:
            self.queue.put((False, "Generating preview"))

        values = self.bounding_box()
        adeviation = max((values[3] - values[0],
                          values[4] - values[1],
                          values[5] - values[2]))

        from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
        BRepMesh_IncrementalMesh(self.shape, 0.05 * adeviation * mesh_quality,
                                 False, 0.5, self.occ_parallel)

        bt = BRep_Tool()

        L = 1 if len(self.faces) == 0 else max(list(self.faces)) + 1
        face_vert_offset = [0] * L

        all_ptx = []
        face_idx = {}
        edge_idx = {}
        vert_idx = {}

        # in order to update value from inner functions. this needs to be
        # object
        offset = Counter()
        num_failedface = Counter()
        num_failededge = Counter()

        solidMap, faceMap, edgeMap, vertMap = self.prep_maps(
            self.shape, return_all=False)

        dprint1("Entity counts: solid/face/edge/vert : ",
                solidMap.Size(), faceMap.Size(), edgeMap.Size(), vertMap.Size())

        solid2isolid = topo2id(self.solids, solidMap)
        face2iface = topo2id(self.faces, faceMap)
        edge2iedge = topo2id(self.edges, edgeMap)
        vert2iverte = topo2id(self.vertices, vertMap)

        face2solid = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            self.shape, TopAbs_FACE, TopAbs_SOLID, face2solid)
        edge2face = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            self.shape, TopAbs_EDGE, TopAbs_FACE, edge2face)
        vertex2edge = TopTools_IndexedDataMapOfShapeListOfShape()
        topexp_MapShapesAndAncestors(
            self.shape, TopAbs_VERTEX, TopAbs_EDGE, vertex2edge)

        def value2coord(value, location):
            if not location.IsIdentity():
                trans = location.Transformation()
                xyz = [v.XYZ() for v in value]
                void = [trans.Transforms(x) for x in xyz]
                ptx = [x.Coord() for x in xyz]
            else:
                ptx = [x.Coord() for x in value]
            return np.vstack(ptx)

        def work_on_face(iface, face):
            face_vert_offset[iface] = offset()

            location = TopLoc_Location()
            facing = (bt.Triangulation(face, location))

            if facing is None:
                num_failedface.increment(1)
                return
            else:
                tab = facing.Nodes()
                tri = facing.Triangles()
                idx = [tri.Value(i).Get()
                       for i in range(1, facing.NbTriangles() + 1)]
                values = [tab.Value(i) for i in range(1, tab.Length() + 1)]
                ptx = value2coord(values, location)

                all_ptx.append(np.vstack(ptx))

                face_idx[iface] = np.vstack(idx) - 1 + offset()
                offset.increment(tab.Length())
                return

        def work_on_edge_on_face(iedge, edge):
            faces = edge2face.FindFromKey(edge)
            topology_iterator = TopTools_ListIteratorOfListOfShape(faces)
            while topology_iterator.More():
                face = topology_iterator.Value()
                topology_iterator.Next()
                location = TopLoc_Location()
                facing = (bt.Triangulation(face, location))
                if facing is not None:
                    break
            else:
                num_failededge.increment(1)
                print('tesselation of edge is missing, iedge=', iedge)
                return

            iface = face2iface[face]
            coffset = face_vert_offset[iface]
            poly = (bt.PolygonOnTriangulation(edge, facing, location))

            if poly is None:
                num_failededge.increment(1)
            else:
                node = poly.Nodes()
                idx = [
                    node.Value(i) +
                    coffset -
                    1 for i in range(
                        1,
                        poly.NbNodes() +
                        1)]
                edge_idx[iedge] = idx

        def work_on_edge(iedge, edge):
            location = TopLoc_Location()
            poly = bt.Polygon3D(edge, location)

            if poly is None:
                work_on_edge_on_face(iedge, edge)
            else:
                nnodes = poly.NbNodes()
                nodes = poly.Nodes()
                values = [nodes.Value(i) for i in range(1, poly.NbNodes() + 1)]
                ptx = value2coord(values, location)

                idx = np.arange(poly.NbNodes())

                all_ptx.append(np.vstack(ptx))
                edge_idx[iedge] = list(idx + offset())
                offset.increment(poly.NbNodes())

        def work_on_vertex(ivert, vertex):
            pnt = bt.Pnt(vertex)
            ptx = [pnt.Coord()]
            idx = [offset()]
            all_ptx.append(ptx)
            vert_idx[ivert] = idx
            offset.increment(1)

        for iface in self.faces:
            work_on_face(iface, self.faces[iface])
        for iedge in self.edges:
            work_on_edge(iedge, self.edges[iedge])
        for ivert in self.vertices:
            work_on_vertex(ivert, self.vertices[ivert])

        def generate_idxmap_from_map(idxmap, parent_imap, child2parents, objs):
            for iobj in objs:
                parents = child2parents.FindFromKey(objs[iobj])
                topology_iterator = TopTools_ListIteratorOfListOfShape(parents)
                while topology_iterator.More():
                    p = topology_iterator.Value()
                    topology_iterator.Next()

                    try:
                        iparent = parent_imap[p]
                    except BaseException:
                        assert False, "Not found"

                    idxmap[iparent].append(iobj)

        # make v, s, l
        v = defaultdict(list)
        s = defaultdict(list)
        l = defaultdict(list)

        generate_idxmap_from_map(v, solid2isolid, face2solid, self.faces)
        generate_idxmap_from_map(s, face2iface, edge2face, self.edges)
        generate_idxmap_from_map(l, edge2iedge, vertex2edge, self.vertices)

        v = dict(v)
        s = dict(s)
        l = dict(l)

        shape = {}
        idx = {}

        # vertex
        keys = list(vert_idx)
        if len(keys) > 0:
            shape['vertex'] = np.vstack([vert_idx[k] for k in keys])
            idx['vertex'] = {'geometrical': np.hstack([k for k in keys]),
                             'physical': np.hstack([0 for k in keys])}
        # edge
        keys = list(edge_idx)
        if len(keys) > 0:
            a = [np.vstack([edge_idx[k][:-1], edge_idx[k][1:]]).transpose()
                 for k in keys]
            shape['line'] = np.vstack(a)
            eidx = np.hstack([[k] * (len(edge_idx[k]) - 1) for k in keys])
            idx['line'] = {'geometrical': eidx,
                           'physical': eidx * 0}

        # face
        keys = list(face_idx)
        if len(keys) > 0:
            shape['triangle'] = np.vstack([face_idx[k] for k in keys])
            eidx = np.hstack([[k] * len(face_idx[k]) for k in keys])
            idx['triangle'] = {'geometrical': eidx,
                               'physical': eidx * 0}

        ptx = np.vstack(all_ptx)
        esize = self.get_esize()

        vcl = self.get_vcl(l, esize)
        geom_msh = ''

        dprint1(
            "number of triangulation fails",
            num_failedface(),
            num_failededge())
        return geom_msh, l, s, v, vcl, esize, ptx, shape, idx

    def generate_brep(self, objs, filename='', trash='', finalize=False):

        if finalize and not self.skip_final_frag:
            if self.logfile is not None:
                self.logfile.write("finalize is on \n")
            if self.queue is not None:
                self.queue.put((False, "finalize is on"))

            self.apply_fragments()
            self.synchronize_topo_list(action='both')
            
        geom_brep = self.make_safe_file(filename, trash, '.brep')
        self.write_brep(geom_brep)

        '''
        do_map_always = False
        if finalize or do_map_always:

            We need to reload it here so that indexing is consistent
            in meshing.

            # We keep highestDimOnly = False, sinse low dim elemtns could be
            # used for embeding (cl control)
            #gmsh.model.occ.importShapes(geom_brep, highestDimOnly=False)
            # gmsh.model.occ.synchronize()
            #self.applyEntityNumberingInfo(map, objs)
        '''
        return geom_brep

    '''
    sequence/preview/brep generator
    '''

    def run_sequence(self, objs, gui_data, start_idx):
        isWP = False

        print("start idx", start_idx)
        if start_idx < 1:
            self.shape = TopoDS_Compound()
            self.builder.MakeCompound(self.shape)
            self.prep_topo_list()

        for gui_name, gui_param, geom_name in self.geom_sequence[start_idx:]:
            if self.logfile is not None:
                self.logfile.write("processing " + gui_name + "\n")
                self.logfile.write(
                    "data " +
                    str(geom_name) +
                    ":" +
                    str(gui_param) +
                    "\n")
            if self.queue is not None:
                self.queue.put((False, "processing " + gui_name))

            if geom_name == "WP_Start":
                tmp = objs.duplicate()
                org_keys = list(objs)

                for x in org_keys:
                    del tmp[x]

                org_objs = objs
                objs = tmp
                isWP = True

            elif geom_name == "WP_End":
                for x in objs:
                    org_objs[x] = objs[x]
                objs = org_objs
                isWP = False

            else:
                try:
                    method = getattr(self, geom_name + '_build_geom')
                    objkeys, newobjs = method(objs, *gui_param)
                    gui_data[gui_name] = (objkeys, newobjs)
                except BaseException:
                    import traceback
                    if self.logfile is not None:
                        self.logfile.write("failed " + traceback.format_exc())
                    assert False, traceback.format_exc()

        #capcheName = "" if isWP else gui_name
        return gui_data, objs


class OCCGeometryGenerator(mp.Process):
    def __init__(self, q, task_q):
        self.q = q
        self.task_q = task_q
        self.mw = None
        mp.Process.__init__(self)

    def run(self):
        while True:
            time.sleep(0.1)
            try:
                task = self.task_q.get(True)
            except EOFError:
                self.result_queue.put((-1, None))
                # self.task_queue.task_done()
                continue

            if task[0] == -1:
                # self.task_queue.task_done()
                break
            if task[0] == 1:
                try:
                    self.generator(*task[1])
                except BaseException:
                    import traceback
                    txt = traceback.format_exc()
                    traceback.print_exc()
                    self.q.put((True, ('fail', txt)))
                    # self.task_queue.task_done()
                    break
        print("exiting prcesss")

    def generator(self, sequence, no_mesh, finalize,
                  filename, start_idx, trash, kwargs):

        kwargs['write_log'] = True
        kwargs['queue'] = self.q
        q = self.q

        if self.mw is None:
            from petram.geom.gmsh_geom_model import GeomObjs
            self.mw = Geometry(**kwargs)
            self.objs = GeomObjs()
            self.gui_data = dict()
        else:
            self.mw.process_kwargs(kwargs)

        q.put((self.mw.logfile.name))

        self.mw.geom_sequence = sequence
        self.mw.run_sequence(self.objs, self.gui_data, start_idx)

        if finalize:
            filename = filename
            brep_file = self.mw.generate_brep(self.objs, filename=filename,
                                              finalize=True)
        else:
            filename = sequence[-1][0]
            brep_file = self.mw.generate_brep(self.objs, filename=filename,
                                              trash=trash, finalize=False)

        if no_mesh:
            q.put((True, (self.gui_data, self.objs, brep_file, None, None)))

        else:
            data = self.mw.generate_preview_mesh(filename, trash)
            # data =  geom_msh, l, s, v,  vcl, esize

            q.put((True, (self.gui_data, self.objs, brep_file, data, None)))
