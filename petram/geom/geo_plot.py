import numpy as np

def plot_geometry(viewer,  ret):
    X, cells, pt_data, cell_data, field_data = ret

    vert = np.squeeze(X[cells['vertex']])
    viewer.plot(vert[:,0],
                vert[:,1],
                vert[:,2], 'o', linewidth = 0)

    lidx = np.unique(cell_data['line']['geometrical'])
    lines = cell_data['line']['geometrical']

    viewer.update(False)
    for i in lidx:
       idx = np.where(lines == i)
       obj = viewer.solid(X[cells['line'][idx]])
       obj.rename('Line_'+str(i))
    viewer.update(True)       
