import numpy as np


def plot_geometry(viewer,  ret):
    viewer.cls()
    X, cells, pt_data, cell_data, field_data = ret
    print cells.keys()
    if 'vertex' in cells:
        vert = np.squeeze(X[cells['vertex']])
        obj= viewer.plot(vert[:,0],
                    vert[:,1],
                    vert[:,2], 'ob',
                    array_idx = cell_data['vertex']['geometrical'],
                    linewidth = 0)
        obj.rename('Points')
        obj._artists[0].set_gl_hl_use_array_idx(True)

    if 'line' in cells:
        obj = viewer.solid(X[cells['line']],
                           array_idx = cell_data['line']['geometrical'],
                           linewidth = 1.5)

        obj.rename('Edges')
        obj._artists[0].set_gl_hl_use_array_idx(True)

    if 'triangle' in cells:
        obj = viewer.solid(X[cells['triangle']],
                           array_idx = cell_data['triangle']['geometrical'],
                           linewidth = 0)

        obj.rename('Faces')
        obj._artists[0].set_gl_hl_use_array_idx(True)
    

