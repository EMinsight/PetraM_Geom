from collections import defaultdict
import numpy as np

gkeys = ['MeshFormat',
         'PhysicalNames',
         'Nodes',
         'Elements',]

from petram.geom.read_gmsh import (gmsh_element_type,
                                   gmsh_element_dim,
                                   gmsh_element_mfemname,
                                   gmsh_element_order,
                                   num_nodes_per_cell,
                                   num_verts_per_cell)
import mfem.ser as mfem

'''
   t = Translator('file.msh')
   t.write('file.mesh')
   t.write_first_order('file.mesh')
   t.generate_nodalspace()

   msh = GmshFile('file.msh')
'''
class Translator():
    def __init__(self, file=None, verbose=False):
        self.msh = None
        if file is not None:
            self.read(file, verbose=verbose)

    def write(self, dst, linear=False, verbose=False):
        mesh = self.generate_mesh(linear=linear, verbose=verbose)
        print('writing file', dst)
        mesh.Print(dst)

    def read(self, file, verbose=False):
        self.msh = GmshFile(file, verbose=verbose)

    def generate_mesh(self, msh=None, linear=False, verbose=False):
        if msh is None:
            msh = self.msh

        sdim = 3
        for i in (3, 2, 1, 0):
            if len(self.msh.elems[i]) != 0:
                elems = self.msh.elems[i]
                belems = self.msh.elems[i-1]
                dim = i
                break

        Nbdrelem = 0
        for k in belems:
            Nbdrelem = len(belems[k])
        Nelem = 0
        for k in elems:
            Nelem = len(elems[k])

        uidx, _r_idx = msh.count_vertex()
        uidx_map = {u: k         for k, u in enumerate(uidx)}
        Nvert = len(uidx)

        mname_dim = {'Hex': 3, 'Tet': 3, 'Wedge': 3,
                     'Triangle': 2, 'Quad': 2,
                     'Segment': 1, 'Vertex': 0}

        if verbose:
            print("Mesh : dim,  Nvert, Nelem, Nbdrelem, sdim : ",
                  dim, Nvert, Nelem, Nbdrelem, sdim)

        mesh = mfem.Mesh(i, Nvert, Nelem, Nbdrelem, sdim)
        for i in uidx:
            pt = list(msh.nodes[i-1])
            mesh.AddVertex(pt)

        for d in msh.elems:
            for kind in d:
                n_nodes = num_nodes_per_cell[kind]
                n_verts = num_verts_per_cell[kind]
                mname = gmsh_element_mfemname[kind]
                if mname_dim[mname] == dim:  # add elemement
                    adder = getattr(mesh, 'Add'+mname)
                elif mname_dim[mname] == dim-1: # add boundary element
                    adder = getattr(mesh, 'AddBdr'+mname)
                else:
                    continue

                for data in d[kind]:
                    verts = data[-n_nodes:(-n_nodes+n_verts)]
                    verts = [uidx_map[v] for v in verts]
                    attr = data[-n_nodes-1]
                    adder(verts, attr)

        mesh.FinalizeMesh()
        if not linear:
            self.add_nodal(mesh, verbose=verbose)
            
        return mesh

    def add_nodal(self, mesh, verbose=False):
        for i in (3, 2, 1, 0):
            if len(self.msh.elems[i]) != 0:
                tmp = [gmsh_element_order[x] for x in self.msh.elems[i]]
                elems = self.msh.elems[i]
                break
        orders = np.unique(tmp)
        if len(orders) != 1:
            assert False, "mesh is mixed element order (not supported)"

        order = orders[0]
        if verbose:
            print("Element oder : ", order)

        sdim = mesh.SpaceDimension()
        fec_type = mfem.H1_FECollection
        fe_coll = fec_type(order, sdim, mfem.BasisType.ClosedUniform)
        nodal_fes = mfem.FiniteElementSpace(mesh, fe_coll, sdim)
        mesh.SetNodalFESpace(nodal_fes)
        mesh._nodal= nodal_fes

        kelem = 0
        for kind in elems:
            kelem = self.move_nodal(mesh, kind, elems[kind], kelem, verbose)

    def move_nodal(self, mesh, kind, elems, kelem, verbose=True):
        msh = self.msh

        nodes = mesh.GetNodes()
        sdim = mesh.SpaceDimension()
        fes = mesh._nodal

        n_nodes = num_nodes_per_cell[kind]
        n_verts = num_verts_per_cell[kind]

        fe = fes.GetFE(kelem)
        ir = fe.GetNodes()

        ipoints = np.array([(ir.IntPoint(i).x, ir.IntPoint(i).y, ir.IntPoint(i).z)
                            for i in range(ir.GetNPoints())])

        # mapping gmsh index to mfem index
        mappers = {'tetra10': [0, 1, 2, 3, 4, 6, 7, 5, 9, 8],
                   'tetra20': [0, 1, 2, 3, 4, 5, 9, 8, 11, 10,
                               6, 7, 15, 14, 13, 12, 19, 18, 17, 16],
                   'tetra35': [0, 1, 2, 3, 4, 5, 6, 12, 11, 10, 15,
                               14, 13, 7, 8, 9, 21, 20, 19, 18, 17,
                               16, 32, 33, 31, 28, 29, 30, 25, 26,
                               27, 22, 23, 24, 34],                   
                   'tetra56': [0, 1, 2, 3, 4, 5, 6, 7, 15, 14, 13, 12,
                               19, 18, 17, 16, 8, 9, 10, 11, 27, 26, 25,
                               24, 23, 22, 21, 20, 47, 50, 48, 49, 51, 46,
                               40, 43, 41, 45, 44, 42, 34, 37, 35, 39, 38,
                               36, 28, 31, 29, 33, 32, 30, 52, 53, 54, 55],
                   'triangle6':[0, 1, 2, 3, 4, 5],
                   'triangle10':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   'triangle15':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                 10, 11, 12, 13, 14],
                   'triangle21':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                 13, 14, 15, 18, 16, 20, 19, 17]}

        if kind in mappers:
            mapper = mappers[kind]
        else:
            print("mapper is not defined for ", kind)
            print("number of integration points", len(ipoints))
            print("integration points", ipoints)
            assert False, "mapper not implemented"

        for data in elems:
            dofs = fes.GetElementDofs(kelem)
            dof_idx = [[fes.DofToVDof(i, d)
                       for d in range(sdim)] for i in dofs]
            verts = data[-n_nodes:]
            for kk in range(len(mapper)):
                kkv = verts[mapper[kk]]-1
                nodes[dof_idx[kk][0]] = msh.nodes[kkv][0]
                nodes[dof_idx[kk][1]] = msh.nodes[kkv][1]
                nodes[dof_idx[kk][2]] = msh.nodes[kkv][2]
            kelem = kelem+1

        return kelem

class GmshFile():
    def __init__(self, file=None, verbose=False):
        if file is not None:
            self.read_file(file, verbose=verbose)

    def load_MeshFormat(self, fid, verbose):
        while True:
            line = fid.readline()
            if line.startswith('$End'):
                break
            self.format = [float(x) for x in line.strip().split(' ')]
        if verbose:
            print("Format :", self.format)

    def load_PhysicalNames(self, fid, verbose):

        self.physicalnames = {}
        line = fid.readline()   # skip first
        while True:
            line = fid.readline()
            if line.startswith('$End'):
                break

            tmp = line.strip().split(' ')
            dimtag = (int(tmp[0]), int(tmp[1]))
            name = tmp[2][1:-1]

            self.physicalnames[dimtag] = name

        #if verbose:
        #    print("Physical Names :", list(self.physicalnames))

    def load_Nodes(self, fid, verbose):
        num_node = int(fid.readline())

        self.nodes = np.zeros((num_node, 3))
        i = 0
        while True:
            line = fid.readline()
            if line.startswith('$End'):
                break
            tmp = line.strip().split(' ')
            self.nodes[i, :] = np.array([float(x) for x in tmp[1:]])
            i = i + 1

        if verbose:
            print("Size of node array: ", self.nodes.shape)

    def load_Elements(self, fid, verbose):
        num_elem = int(fid.readline())

        self.elems = [defaultdict(list),
                      defaultdict(list),
                      defaultdict(list),
                      defaultdict(list),]

        i = 0
        while True:
            line = fid.readline()
            if line.startswith('$End'):
                break

            tmp = [int(x) for x in line.strip().split(' ')]
            kind = gmsh_element_type[tmp[1]]
            data = tmp[2:]  # data is number-of-phys, phys, phys, nodes
            dim = gmsh_element_dim[tmp[1]]

            self.elems[dim][kind].append(data)
            i = i + 1

        if verbose:
            print("number of element read", i, "out of", num_elem)

        for k in range(4):
            self.elems[k] = dict(self.elems[k])

    def read_file(self, file, verbose=False):
        fid = open(file)

        while True:
            l = fid.readline().strip()
            if l == '': break
            if l.startswith("$"):
                sc = l[1:]
                if verbose:
                    print("Reading :", sc)
                m = getattr(self, 'load_'+sc)
                m(fid, verbose)
        if verbose:
            print("Found ... ", [list(d) for d in self.elems])
            print("Done.")
        fid.close()

    def count_vertex(self):
        idx = []
        for d in self.elems:
            for kind in d:
                n_nodes = num_nodes_per_cell[kind]
                n_verts = num_verts_per_cell[kind]
                for data in d[kind]:
                    verts = data[-n_nodes:(-n_nodes+n_verts)]
                    idx.append(verts)
        idx = np.hstack(idx)

        u_idx, r_idx =  np.unique(idx, return_inverse=True)
        return u_idx, r_idx 
