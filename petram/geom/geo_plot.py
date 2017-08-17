import numpy as np

def expand_vertex_data(X, vertex_idx, element_id):
    '''
    expand index data using element_id, so that
   
       surface edges will have duplicate vertex,
       this way,,
         1) the normal vector on the edge become discutinous
         2) the elment index on the surface becomes constant
    '''
    k = 0
    verts = []
    iele = []
    iarr = []
        
    nel = vertex_idx.shape[-1]
    for kk in np.unique(element_id):
        idx = np.where(element_id == kk)[0]

        iverts = vertex_idx[idx].flatten()
        iv, idx = np.unique(iverts, return_inverse = True)
        verts.append(X[iv])
        iele.append(idx.reshape(-1, nel)+k)
        k = k + len(iv)
        iarr.append(np.zeros(len(iv))+kk)

    array_idx = np.hstack(iarr)
    elem_idx = np.vstack(iele)
    verts = np.vstack(verts)

    return verts, elem_idx, array_idx

def plot_geometry(viewer,  ret):
    viewer.cls()
    viewer.get_axes()._artists[0].gl_hl_setcolor([1,0,0,])
    
    X, cells, pt_data, cell_data, field_data = ret
    print cells.keys()
    if 'vertex' in cells:
        vert = np.squeeze(X[cells['vertex']])
        obj= viewer.plot(vert[:,0],
                    vert[:,1],
                    vert[:,2], 'ob',
                    array_idx = cell_data['vertex']['geometrical'],
                    linewidth = 0)
        obj.rename('points')
        obj._artists[0].set_gl_hl_use_array_idx(True)
 
    if 'line' in cells:
        verts, elem_idx, array_idx = expand_vertex_data(X, cells['line'],
                                       cell_data['line']['geometrical'])
        
        obj = viewer.solid(verts, elem_idx, 
                           array_idx = array_idx,
                           linewidth = 1.5)

        obj.rename('edges')
        obj._artists[0].set_gl_hl_use_array_idx(True)

    if 'triangle' in cells:
        verts, elem_idx, array_idx = expand_vertex_data(X, cells['triangle'],
                                       cell_data['triangle']['geometrical'])


        #print verts.shape, elem_idx.shape, array_idx.shape
        obj = viewer.solid(verts, elem_idx,
                           array_idx = array_idx,
                           facecolor = 'cyan',
                           linewidth = 0)

        obj.rename('faces')
        obj._artists[0].set_gl_hl_use_array_idx(True)

    viewer.set_sel_mode(viewer.get_sel_mode())

