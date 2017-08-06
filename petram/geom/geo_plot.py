import numpy as np

def plot_geometry(viewer,  ret):
    viewer.cls()
    X, cells, pt_data, cell_data, field_data = ret
    if not 'vertex' in cells: return
    vert = np.squeeze(X[cells['vertex']])
    obj= viewer.plot(vert[:,0],
                vert[:,1],
                vert[:,2], 'ob',
                array_idx = cell_data['vertex']['geometrical'],
                linewidth = 0)
    obj.rename('Points')
    obj._artists[0].set_gl_hl_use_array_idx(True)
    
    lidx = np.unique(cell_data['line']['geometrical'])
    lines = cell_data['line']['geometrical']

    obj = viewer.solid(X[cells['line']],
                       array_idx = cell_data['line']['geometrical'],
                       linewidth = 1.5)

    obj.rename('Edges')
    obj._artists[0].set_gl_hl_use_array_idx(True)

