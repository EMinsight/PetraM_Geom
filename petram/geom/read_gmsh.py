import numpy as np

gmsh_element_type = {
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
        15: 'vertex',
        16: 'quad8',
        17: 'hexahedron20',
        18: 'prism15',
        19: 'pyramid13',
        20: 'itriangle9',   # incomplete-triangle
        21: 'triangle10', 
        22: 'itriangle12',  # incomplete-triangle
        23: 'triangle15',
        24: 'itriangle15',  # incomplete-triangle
        25: 'triangle21',
        26: 'line4',    
        27: 'line5',    
        28: 'line6',    
        29: 'tetra20',
        30: 'tetra35',
        31: 'tetra56',
        92: 'hexahedron64',
        93: 'hexahedron125',
        }

gmsh_element_dim = {
        1:  1,
        2:  2,
        3:  2,
        4:  3,
        5:  3,
        6:  3,
        7:  3,
        8:  1,
        9:  2,
        10: 2,
        11: 3,
        12: 3,
        13: 3,
        14: 3,
        15: 0,
        16: 2,
        17: 3,
        18: 3,
        19: 3,
        20: 2,
        21: 2,
        22: 2,
        23: 2,
        24: 2,
        25: 2,
        26: 1,
        27: 1,
        28: 1,
        29: 3,
        30: 3,
        31: 3,
        92: 3,
        93: 3,
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

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('ReadGMSH')

#dimtags =  gmsh.model.getEntities()
def read_loops(geom, dimtags=None):
    model = geom.model
    
    model.occ.synchronize()
    v = {}
    s = {}
    l = {}

    dimtag3 =  model.getEntities(3) if dimtags is None else [x for x in dimtags if x[0]==3]
    for dim, tag in dimtag3:
        v[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]
    dimtag2 =  model.getEntities(2) if dimtags is None else [x for x in dimtags if x[0]==2] + model.getBoundary(dimtag3, combined=False, oriented=False)
    dimtag2 = list(set(dimtag2))
    
    for dim, tag in dimtag2:
        s[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]

    dimtag1 =  model.getEntities(1) if dimtags is None else [x for x in dimtags if x[0]==1] + model.getBoundary(dimtag2, combined=False, oriented=False)
    dimtag1 = list(set(dimtag1))
    for dim, tag in dimtag1:
        l[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False, recursive=True)]
        if len(l[tag]) == 0:
             print('line :' + str(tag) + ' has no boundary (loop)')
             node, coord, pc  = model.mesh.getNodes(dim=1, tag=tag, includeBoundary=True)
             l[tag].append(node[0])
        
    return l, s, v

def read_loops_do_meshloop(geom, dimtags=None):
    model = geom.model
    model.occ.synchronize()
    ent = model.getEntities()    
    model.setVisibility(ent, False)
    
    v = {}
    s = {}
    l = {}

    dimtag3 =  model.getEntities(3) if dimtags is None else [x for x in dimtags if x[0]==3]
    for dim, tag in dimtag3:
        v[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]
    dimtag2 =  model.getEntities(2) if dimtags is None else [x for x in dimtags if x[0]==2] + model.getBoundary(dimtag3, combined=False, oriented=False)
    dimtag2 = list(set(dimtag2))
    
    for dim, tag in dimtag2:
        s[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False)]

    dimtag1 =  model.getEntities(1) if dimtags is None else [x for x in dimtags if x[0]==1] + model.getBoundary(dimtag2, combined=False, oriented=False)
    dimtag1 = list(set(dimtag1))
    for dim, tag in dimtag1:
        l[tag] = [y for x, y in model.getBoundary([(dim, tag)],
                                                       oriented=False, recursive=True)]

        if len(l[tag]) == 0:
             print('line :' + str(tag) + ' has no boundary (loop)')
             model.setVisibility(((dim,tag),), True, recursive = True)
             model.mesh.generate(1)
        
             node, coord, pc  = model.mesh.getNodes(dim=1, tag=tag, includeBoundary=True)
             l[tag].append(node[0])
        
    return l, s, v



def read_loops2(geom, dimtags=None):
    '''
    read vertex coordinats and loops together. 
    before calling this, do 
        self.hide_all()        
        gmsh.model.mesh.generate(1)
    '''
    l, s, v = read_loops(geom, dimtags=dimtags)
    model = geom.model

    nidx, coord, pcoord = geom.model.mesh.getNodes(dim=0)
    ptx = np.array(coord).reshape(-1, 3)
    
    if dimtags is not None:
        ents = model.getBoundary(dimtags, combined=False, recursive=True)
        tags = list(set([tag for didm, tag in ents]))
        p = {t:int(geom.model.mesh.getNodes(dim=0, tag=t)[0][0]-1) for t in tags}
    else:
        tags = [tag for dim, tag in geom.model.getEntities(0)]
        p = {t: int(nidx[k]-1)  for k, t in enumerate(tags)}        

    return ptx, p, l, s, v

def read_loops3(geom, dimtags=None):
    ptx, p, l, s, v = read_loops2(geom, dimtags=dimtags)

    mid_points = {}
    for k in l:
        p1 = ptx[p[l[k][0]]]
        p2 = ptx[p[l[k][1]]]
        param1 = geom.model.getParametrization(1, k, p1)
        param2 = geom.model.getParametrization(1, k, p2)       
        value = geom.model.getValue(1, k, (param1 + param2)/2.0)
        mid_points[k] = value
        #print(param1, param2, value)

    return ptx, p, l, s, v, mid_points

def read_pts_groups(geom, finished_lines=None, 
                          finished_faces=None):

    # w/o using this it becomes so slow...
    finished_lines=None
    finished_faces=None

    dprint1("finished faces", finished_faces)
    dprint1("finished lines", finished_lines)    
    #finished_faces is None: finished_faces = []     
    model = geom.model
    
    node_id, node_coords, parametric_coods =  model.mesh.getNodes()

    if len(node_coords) == 0:
        return np.array([]).reshape((-1,3)), {}, {}
    points = np.array(node_coords).reshape(-1, 3)

    node2idx = np.zeros(int(max(node_id)+1), dtype=int)

    for k, id in enumerate(node_id): node2idx[id] = k

    # cells is element_type -> node_id 
    cells = {}
    cell_data = {}
    el2idx = {}
    for ndim in range(3):
        dimtags =  model.getEntities(dim=ndim)
        
        if (ndim == 1 and finished_lines is not None and len(finished_lines) != 0
            or
            ndim == 2 and finished_faces is not None and len(finished_faces) != 0):
            finished = finished_faces if ndim==2 else finished_lines

            dimtags = [dt for dt in dimtags if dt[1] in finished]
            
            xxx = [model.mesh.getElements(ndim, l[1]) for l in dimtags]

            if len(xxx) == 0: continue
            if isinstance(xxx[0][0], np.ndarray):
                # newer gmsh comes here
                elementTypes = np.hstack([x[0] for x in xxx])
                elementTags = sum([x[1] for x in xxx],[])
                nodeTags = sum([x[2] for x in xxx],[])
                elementTags = [list(x) for x in elementTags]
                nodeTags = [list(x) for x in nodeTags]                
            else:

                elementTypes = sum([x[0] for x in xxx], [])
                elementTags =sum([x[1] for x in xxx], [])
                nodeTags = sum([x[2] for x in xxx], [])
                
            tmp = {}                
            dd1 = {};
            dd2 = {}
            for k, el_type in enumerate(elementTypes):
                if not el_type in dd1: dd1[el_type] = []
                if not el_type in dd2: dd2[el_type] = []
                dd1[el_type] = dd1[el_type]+ elementTags[k]
                dd2[el_type] = dd2[el_type]+ nodeTags[k]
            elementTypes = list(dd1)
            elementTags = [dd1[k] for k in elementTypes]
            nodeTags = [dd2[k] for k in elementTypes]
        else:
            elementTypes, elementTags, nodeTags = model.mesh.getElements(ndim)

        for k, el_type in enumerate(elementTypes):
            el_type_name = gmsh_element_type[el_type]
            data = np.array([node2idx[tag] for tag in nodeTags[k]], dtype=int)
            data = data.reshape(-1, num_nodes_per_cell[el_type_name])
            cells[el_type_name] = data
            tmp = np.zeros(int(max(elementTags[k])+1), dtype=int)
            for kk, id in enumerate(elementTags[k]): tmp[id] = kk
            el2idx[el_type_name] = tmp
            cell_data[el_type_name] = {'geometrical':
                                       np.zeros(len(elementTags[k]), dtype=int),
                                       'physical':
                                       np.zeros(len(elementTags[k]), dtype=int)}
        '''
        if (ndim == 1 and finished_lines is not None
            or
            ndim == 2 and finished_faces is not None):
            finished = finished_faces if ndim==2 else finished_lines
            #print("here we are removing unfinished lines",  dimtags, finished_lines)
            dimtags = [dt for dt in dimtags if dt[1] in finished]
        '''
        for dim, tag in dimtags:
            elType2, elTag2, nodeTag2 = model.mesh.getElements(dim=dim,
                                                               tag=tag)
            for k, el_type in enumerate(elType2):                       
                el_type_name = gmsh_element_type[el_type]
                for elTag in elTag2[k]:
                   if not el_type_name in el2idx: continue
                   idx = el2idx[el_type_name][elTag]
                   cell_data[el_type_name]['geometrical'][idx] = tag
        '''
        dimtags = model.getPhysicalGroups(dim=ndim)
        print("physical", dimtags)
        for dim, ptag in dimtags:
            etags = model.getEntitiesForPhysicalGroup(dim=dim, tag=ptag)
            for etag in etags:
                elType2, elTag2, nodeTag2 = model.mesh.getElements(dim=dim,
                                                                   tag=etag)
                for k, el_type in enumerate(elType2):                       
                    el_type_name = gmsh_element_type[el_type]
                    for elTag in elTag2[k]:
                        if not el_type_name in el2idx: continue                        
                        idx = el2idx[el_type_name][elTag]
                        cell_data[el_type_name]['physical'][idx] = ptag
        '''
    return points, cells, cell_data
