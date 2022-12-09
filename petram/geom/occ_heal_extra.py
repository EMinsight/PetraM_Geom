'''
 OCC heal shape extra
'''
from petram.geom.occ_cbook import *
from scipy.spatial import distance_matrix


def split_face_extra(shape, limit=0.1, tolerance=1e-6):
    '''
    x------------------x-----------------x
    |                                    |
    x------------------x-----------------x
    '''
    bt = BRep_Tool()

    vertices = [p for p in iter_shape_once(shape, 'vertex')]

    ptx = []
    for v in vertices:
        pnt = bt.Pnt(v)
        p = np.array((pnt.X(), pnt.Y(), pnt.Z(),))
        ptx.append(p)

    ptx = np.vstack(ptx)
    print(ptx.shape)
    md = distance_matrix(ptx, ptx, p=2)

    short_pair = []
    for i in range(len(ptx)):
       for j in range(len(ptx)):
           if j > i: break
           if md[i, j] < limit:
               if np.abs(i-j) == 0: continue
               if np.abs(i-j) == 1: continue
               if i == 0 and j == len(ptx)-1: continue
               short_pair.append((i, j))
    print("short distance pair", short_pair)

    assert len(
        short_pair) < 2, "more than two short pairs are found (try largeer limit?)"

    p1 = min(short_pair[0])
    p2 = max(short_pair[0])

    # make a dict to follow edges
    edges1 = {i: i+1 for i in range(len(ptx))}
    edges1[len(ptx)-1] = 0
    print(edges1)
    loop1 = [p1]
    while True:
        pp = edges1[loop1[-1]]
        loop1.append(pp)
        if pp == p2:
            # loop1.append(p1)
            break

    loop2 = [p2]
    while True:
        pp = edges1[loop2[-1]]
        loop2.append(pp)
        if pp == p1:
            # loop2.append(p2)
            break

    print("loop1", loop1)
    print("loop2", loop2)

    # creat a new edge
    edgeMaker = BRepBuilderAPI_MakeEdge(vertices[p1], vertices[p2])
    edgeMaker.Build()
    if not edgeMaker.IsDone():
        assert False, "Can not make line"
    new_edge = edgeMaker.Edge()

    edge_connection = {}
    edges = [p for p in iter_shape_once(shape, 'edge')]
    for e in edges:
        mapper = get_mapper(e, 'vertex')
        idx = np.where([mapper.Contains(v) for v in vertices])[0]
        edge_connection[tuple(idx)] = e

    # make a list of edges to create surface
    edges1 = [edge_connection[(loop1[i], loop1[i+1])]
                               for i in range(len(loop1)-1)]
    edges1.append(new_edges1)
    edges2 = [edge_connection[(loop2[i], loop2[i+1])]
                               for i in range(len(loop2)-1)]
    edges2.append(new_edges1)

    def make_wire(edges):
        wireMaker = BRepBuilderAPI_MakeWire()
        for e in edges:
            wireMaker.Add(e)
            wireMaker.Build()
         return wireMaker.Wire()
         retu

    wire1 = make_wire(edges1)
    wire2 = make_wire(edges2)

    def make_filling(wire):
        f = BRepOffsetAPI_MakeFilling()
        # make wire constraints
        ex1 = BRepTools_WireExplorer(wire)
        while ex1.More():
            edge = topods_Edge(ex1.Current())
            f.Add(edge, GeomAbs_C0)
            ex1.Next()
        f.Build()

        if not f.IsDone():
            assert False, "Cannot make filling"
            
        face = f.Shape()
        s = bt.Surface(face)

        faceMaker = BRepBuilderAPI_MakeFace(s, wire)
        result = faceMaker.Face()

        fix = ShapeFix_Face(result)
        fix.SetPrecision(self.occ_geom_tolerance)
        fix.Perform()
        fix.FixOrientation()
        result = fix.Face()
    
    
    
        from OCC.Core.GeomAbs import GeomAbs_C0
        from OCC.Core.BRepTools import BRepTools_WireExplorer


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
    

        
    '''
