'''

   Direct OCC Brep data reader
   
   This routine can be tested using
      python -c <this file>

'''
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.TopExp import TopExp_Explorer, topexp_MapShapes, topexp_MapShapesAndAncestors
from OCC.Core.TopoDS import TopoDS_Compound, topods_Face, TopoDS_Shape, topods_Edge, topods_Vertex, topods_Solid, topods_Wire, topods_Shell
from OCC.Core.TopAbs import TopAbs_SOLID, TopAbs_WIRE, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SHELL
from OCC.Core.TopTools import TopTools_IndexedMapOfShape
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepTools import breptools_Read

import numpy as np
import os

class TopoShapeSet(list):
    def __init__(self, mapping=None):
        self.mapping = mapping
        if mapping is not None:
            self.check = np.array([-1]*mapping.Size())
            
    def add_shape(self, x):
        if self.mapping is not None:
            i = self.mapping.FindIndex(x)-1
            if self.check[i] != -1:
                return False, self.check[i]+1
            else:
                self.check[i]=len(self)
                self.append(x)
                return True, self.check[i]+1
        else:
            for kk, i in enumerate(self):
                if i.IsSame(x): return False, kk+1
            self.append(x)
            return True, len(self)

    def check_shape(self, x):
        if self.mapping is not None:
            i = self.mapping.FindIndex(x)-1
            if self.check[i] != -1:
                return False
            else:
                return True
        else:
            for kk, i in enumerate(self):
                if i.IsSame(x): return False
            return True
        
class Counter(object):
    def __init__(self):
        self.value = 0
        
    def increment(self,x):
        self.value = self.value + x

    def __call__(self):
        return self.value
    
def read_file(filename, verbose = False):
    if not os.path.exists(filename):
        assert False, "file does not exists"
    #
    # Create the shape
    #
    shape = TopoDS_Shape()
    builder = BRep_Builder()

    success = breptools_Read(shape, filename, builder)

    if not success:
        assert False, "Falied to read file"

    #
    # Mesh the shape
    #
    BRepMesh_IncrementalMesh(shape, 0.5)

    comp = TopoDS_Compound()
    builder.MakeCompound(comp)

    bt = BRep_Tool()    

    all_ptx = []
    face_idx = {}
    edge_idx = {}
    vert_idx = {}
    
    offset = Counter() ## in order to update value from inner functions. this needs to be object
    num_failedface=Counter()
    num_failededge=Counter()

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

    if verbose:
        print("Number of solid/shell/face/wire/edge", solidMap.Size(),
              shellMap.Size(), faceMap.Size(), wireMap.Size(), edgeMap.Size())

    ushells = TopoShapeSet(mapping=shellMap)          
    ufaces = TopoShapeSet(mapping=faceMap)          
    uwires = TopoShapeSet(mapping=wireMap)          
    uedges = TopoShapeSet(mapping=edgeMap)          
    uvertices = TopoShapeSet(mapping=vertMap)
    
    def work_on_shell(shell):
        flag, ishell= ushells.add_shape(shell)
        if not flag: return
        
    def work_on_face(face):
        coffset = offset()
        
        flag, iface = ufaces.add_shape(face)
        if not flag: return None, coffset

        location = TopLoc_Location()
        facing = (bt.Triangulation(face, location))
        trans = location.Transformation()        
        if facing is None:
            num_failedface.increment(1)
            return None, coffset
        else:
            tab = facing.Nodes()
            tri = facing.Triangles()
            idx = [tri.Value(i).Get() for i in range(1, facing.NbTriangles()+1)]
            ptx = []
            for i in range(1, tab.Length()+1):
                p = tab.Value(i).XYZ()
                trans.Transforms(p)
                ptx.append(p.Coord())

            all_ptx.append(np.vstack(ptx))
            
            face_idx[iface] = np.vstack(idx)-1 + offset()
            offset.increment(tab.Length())
        return facing, coffset
    
    def work_on_edge_on_face(edge, facing, coffset):
        flag = uedges.check_shape(edge)
        if not flag: return
        if facing is None: return
        
        location = TopLoc_Location()
        poly = (bt.PolygonOnTriangulation(edge, facing, location))
        if poly is None:
            num_failededge.increment(1)            
            idx = []
        else:
            node = poly.Nodes()
            idx = [node.Value(i)+coffset-1  for i in range(1, poly.NbNodes()+1)]
            flag, iedge = uedges.add_shape(edge)
            edge_idx[iedge] = idx

    def work_on_wire(wire):
        flag, iwire = uwires.add_shape(wire)
        if not flag: return
        
    def work_on_edge(edge):
        flag = uedges.check_shape(edge)
        if not flag: return

        location = TopLoc_Location()
        poly=bt.Polygon3D(edge, location)
        trans = location.Transformation()        
        if poly is None:
            pass
        else:
            nnodes = poly.NbNodes()
            nodes = poly.Nodes()
            ptx = []
            for i in range(1, poly.NbNodes()+1):
                p = nodes.Value(i).XYZ()
                trans.Transforms(p)
                ptx.append(p.Coord())
            
            idx = np.arange(poly.NbNodes())
           
            all_ptx.append(np.vstack(ptx))
            flag, iedge = uedges.add_shape(edge)            
            edge_idx[iedge] = np.vstack(idx) + offset()
            offset.increment(poly.NbNodes())
        
    def work_on_vertex(vertex):
        flag, ivert = uvertices.add_shape(vertex)
        if not flag: return
        
        pnt =bt.Pnt(vertex)
        ptx = [pnt.Coord()]
        idx = [offset()]
        all_ptx.append(ptx)
        vert_idx[ivert] = idx
        offset.increment(1)


    ex1 = TopExp_Explorer(shape, TopAbs_SOLID)
    while ex1.More():
        solid = topods_Solid(ex1.Current())
        ex2 = TopExp_Explorer(solid, TopAbs_SHELL)
        while ex2.More():
            shell = topods_Shell(ex2.Current())
            work_on_shell(shell)
            ex3 = TopExp_Explorer(shell, TopAbs_FACE)
            while ex3.More():    
                face = topods_Face(ex3.Current())
                facing, coffset = work_on_face(face)
                ex4 = TopExp_Explorer(face, TopAbs_WIRE)
                while ex4.More():
                    wire = topods_Wire(ex4.Current())
                    work_on_wire(wire)                                        
                    ex5 = TopExp_Explorer(wire, TopAbs_EDGE)
                    while ex5.More():
                        edge = topods_Edge(ex5.Current())
                        work_on_edge_on_face(edge, facing, coffset)
                        ex6 = TopExp_Explorer(edge, TopAbs_VERTEX)
                        while ex6.More():
                            vertex = topods_Vertex(ex6.Current())
                            work_on_vertex(vertex)
                            ex6.Next()
                        ex5.Next()
                    ex4.Next()
                ex3.Next()
            ex2.Next()
        ex1.Next()
        

    ex2.Init(shape, TopAbs_SHELL, TopAbs_SOLID)
    while ex2.More():
        shell = topods_Shell(ex2.Current())
        flag = ushells.check_shape(shell)        
        if flag:
            work_on_shell(shell)            
            ex3 = TopExp_Explorer(shell, TopAbs_FACE)
            while ex3.More():    
                face = topods_Face(ex3.Current())
                facing, coffset = work_on_face(face)
                ex4 = TopExp_Explorer(face, TopAbs_WIRE)
                while ex4.More():
                    wire = topods_Wire(ex4.Current())
                    work_on_wire(wire)                    
                    ex5 = TopExp_Explorer(wire, TopAbs_EDGE)
                    while ex5.More():
                        edge = topods_Edge(ex5.Current())
                        work_on_edge_on_face(edge, facing, coffset)
                        ex6 = TopExp_Explorer(edge, TopAbs_VERTEX)
                        while ex6.More():
                            vertex = topods_Vertex(ex6.Current())
                            work_on_vertex(vertex)
                            ex6.Next()
                        ex5.Next()
                    ex4.Next()
                ex3.Next()
        ex2.Next()


    ex3 = TopExp_Explorer(shape, TopAbs_FACE, TopAbs_SHELL)
    while ex3.More():
        face = topods_Face(ex3.Current())
        flag = ufaces.check_shape(face)                
        if flag:
            facing, coffset = work_on_face(face)
            ex4 = TopExp_Explorer(face, TopAbs_WIRE)
            while ex4.More():
                wire = topods_Wire(ex4.Current())
                work_on_wire(wire)
                ex5 = TopExp_Explorer(wire, TopAbs_EDGE)
                while ex5.More():
                    edge = topods_Edge(ex5.Current())
                    work_on_edge_on_face(edge, facing, coffset)
                    ex6 = TopExp_Explorer(edge, TopAbs_VERTEX)
                    while ex6.More():
                        vertex = topods_Vertex(ex6.Current())
                        work_on_vertex(vertex)
                        ex6.Next()
                    ex5.Next()
                ex4.Next()
        ex3.Next()

    ex4 = TopExp_Explorer(shape, TopAbs_WIRE, TopAbs_FACE)
    while ex4.More():
        wire = topods_Wire(ex4.Current())
        flag = uwires.check_shape(wire)                        
        if flag:
            work_on_wire(wire)
            ex5 = TopExp_Explorer(wire, TopAbs_EDGE)
            while ex5.More():
                edge = topods_Edge(ex5.Current())
                work_on_edge(edge)
                ex6 = TopExp_Explorer(edge, TopAbs_VERTEX)
                while ex6.More():
                    vertex = topods_Vertex(ex6.Current())
                    work_on_vertex(vertex)
                    ex6.Next()
                ex5.Next()
        ex4.Next()
      

    ex5 = TopExp_Explorer(shape, TopAbs_EDGE, TopAbs_WIRE)
    while ex5.More():
        edge = topods_Edge(ex5.Current())
        flag = uedges.check_shape(edge)
        if flag:
            work_on_edge(edge)
            ex6 = TopExp_Explorer(edge, TopAbs_VERTEX)
            while ex6.More():
                vertex = topods_Vertex(ex6.Current())
                work_on_vertex(vertex)
                ex6.Next()
        ex5.Next()

    ex6 = TopExp_Explorer(shape, TopAbs_VERTEX, TopAbs_EDGE)
    while ex6.More():
        vertex = topods_Vertex(ex6.Current())
        flag  = uvertices.check_shape(vertex)
        if flag:
             work_on_vertex(vertex)
        ex6.Next()

    ptx = np.vstack(all_ptx)
    print(edge_idx)
    return ptx, face_idx, edge_idx, vert_idx, num_failedface, num_failededge


def format_brepdata(dataset, mappings):
    '''
    format data read from brep for plotting
    '''
    ptx = dataset[0]
    face_idx = dataset[1]
    edge_idx = dataset[2]
    vert_idx = dataset[3]
    
    dataset = ptx, face_idx, edge_idx, vert_idx, mappings
    vmap, fmap, emap, pmap = mappings

    shape = {}
    idx = {}

    # vertex
    keys = list(vert_idx)
    if len(keys) > 0:
        shape['vertex'] = np.vstack([vert_idx[k] for k in keys])
        idx['vertex']    = {'geometrical': np.hstack([pmap[k] for k in keys]),
                        'physical': np.hstack([0 for k in keys])}

    # edge
    keys = list(edge_idx)
    if len(keys) > 0:    
        shape['line'] = np.vstack([np.vstack([edge_idx[k][:-1], edge_idx[k][1:]]).transpose()
                                for k in keys])
        eidx = np.hstack([[emap[k]]*(len(edge_idx[k])-1) for k in keys])
        idx['line']    = {'geometrical': eidx,
                        'physical': eidx*0}

    # face
    keys = list(face_idx)
    if len(keys) > 0:
        shape['triangle'] = np.vstack([face_idx[k] for k in keys])
        eidx = np.hstack([[fmap[k]]*len(face_idx[k]) for k in keys])
        idx['triangle']    = {'geometrical': eidx,
                              'physical': eidx*0}

        

    return ptx, shape, idx


if __name__ == '__main__':
    import cProfile, pstats, io
    from pstats import SortKey
    pr = cProfile.Profile()
    pr.enable()
    
    #filename = '/Users/shiraiwa/Desktop/check_numbering.brep'
    #filename = '/Users/shiraiwa/Desktop/test.brep'
    #filename = '/Users/shiraiwa/Desktop/ILA.brep'
    #filename = '/Users/shiraiwa/Desktop/line.brep'
    #filename = '/Users/shiraiwa/Desktop/extrude.brep'
    #filename = '/Users/shiraiwa/Desktop/box.brep'
    filename = '/Users/shiraiwa/Desktop/NSTX_1s_cmplx_shield.brep'
    
    ptx, face_idx, edge_idx, vert_idx, num_failedface, num_failededge = read_file(filename, verbose=True)
    print("number of triangulation fails", num_failedface(), num_failededge())
    pr.dump_stats("test4.prof")


                     
    

