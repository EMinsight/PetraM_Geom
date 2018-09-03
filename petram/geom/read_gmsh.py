import numpy as np

gmsh_element_type = {
        15: 'vertex',
        1: 'line',
        2: 'triangle',
        3: 'quad',
        4: 'tetra',
        5: 'hexahedron',
        6: 'wedge',
        7: 'pyramid',
        8: 'line3',
        9: 'triangle6',
        10: 'quad9',
        11: 'tetra10',
        12: 'hexahedron27',
        13: 'prism18',
        14: 'pyramid14',
        26: 'line4',
        36: 'quad16',
        }
num_nodes_per_cell = {
    'vertex': 1,
    'line': 2,
    'triangle': 3,
    'quad': 4,
    'tetra': 4,
    'hexahedron': 8,
    'wedge': 6,
    'pyramid': 5,
    #
    'line3': 3,
    'triangle6': 6,
    'quad9': 9,
    'tetra10': 10,
    'hexahedron27': 27,
    'prism18': 18,
    'pyramid14': 14,
    'line4': 4,
    'quad16': 16,
    }

#dimtags =  gmsh.model.getEntities()
def read_loops(geom):
    model = geom.model
    
    model.occ.synchronize()
    v = {}
    s = {}
    l = {}
    
    dimtags =  model.getEntities(3)
    for dim, tag in dimtags:
        v[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]
    dimtags =  model.getEntities(2)
    for dim, tag in dimtags:
        s[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]
    dimtags =  model.getEntities(1)
    for dim, tag in dimtags:
        l[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]
    return l, s, v

def read_pts_groups(geom):

    model = geom.model
    
    node_id, node_coords, parametric_coods =  model.mesh.getNodes()
    if len(node_coords) == 0:
        return np.array([]).reshape((-1,3)), {}, {}
    points = np.array(node_coords).reshape(-1, 3)

    node2idx = np.zeros(node_id[-1]+1, dtype=int)
    for k, id in enumerate(node_id): node2idx[id] = k

    # cells is element_type -> node_id 
    cells = {}
    cell_data = {}
    el2idx = {}
    for ndim in range(3):
        elementTypes, elementTags, nodeTags = model.mesh.getElements(ndim)
        for k, el_type in enumerate(elementTypes):
            el_type_name = gmsh_element_type[el_type]
            data = np.array([node2idx[tag] for tag in nodeTags[k]], dtype=int)
            data = data.reshape(-1, num_nodes_per_cell[el_type_name])
            cells[el_type_name] = data

            elementTags
            tmp = np.zeros(max(elementTags[k])+1, dtype=int)
            for kk, id in enumerate(elementTags[k]): tmp[id] = kk
            el2idx[el_type_name] = tmp
            cell_data[el_type_name] = {'geometrical':
                                       np.zeros(len(elementTags[k]), dtype=int),
                                       'physical':
                                       np.zeros(len(elementTags[k]), dtype=int)}

        dimtags =  model.getEntities(dim=ndim)
        for dim, tag in dimtags:
            elType2, elTag2, nodeTag2 = model.mesh.getElements(dim=dim,
                                                                    tag=tag)
            for k, el_type in enumerate(elType2):                       
                el_type_name = gmsh_element_type[el_type]
                for elTag in elTag2[k]:
                   idx = el2idx[el_type_name][elTag]
                   cell_data[el_type_name]['geometrical'][idx] = tag

        dimtags = model.getPhysicalGroups(dim=dim)
        for dim, ptag in dimtags:
            etags = model.getEntitiesForPhysicalGroup(dim=dim, tag=ptag)
            for etag in etags:
                elType2, elTag2, nodeTag2 = model.mesh.getElements(dim=dim,
                                                                        tag=etag)
                for k, el_type in enumerate(elType2):                       
                    el_type_name = gmsh_element_type[el_type]
                    for elTag in elTag2[k]:
                        idx = el2idx[el_type_name][elTag]
                        cell_data[el_type_name]['physical'][idx] = ptag


    return points, cells, cell_data
