import numpy as np

from petram.utils import get_pkg_datafile
import petram.geom

fdot = get_pkg_datafile(petram.pi, 'icon',  'dot.png')
fedge = get_pkg_datafile(petram.pi, 'icon', 'line.png')
fface = get_pkg_datafile(petram.pi, 'icon', 'face.png')
fdom = get_pkg_datafile(petram.pi, 'icon', 'domain.png')
show = get_pkg_datafile(petram.pi, 'icon', 'show.png')
hide = get_pkg_datafile(petram.pi, 'icon', 'hide.png')

from petram.pi.sel_buttons import _select_x

def select_dot(evt):
    _select_x(evt, 'point', 'point')
    
def select_edge(evt):
    _select_x(evt, 'edge', 'edge')
    
def select_face(evt):
    _select_x(evt, 'face', 'face')
    
def select_volume(evt):
    _select_x(evt, 'volume', 'face')    

def show_all(evt):
    viewer = evt.GetEventObject().GetTopLevelParent()
    mode = viewer._sel_mode

    ax = viewer.get_axes()
    if mode == 'volume':
        ax.face.hide_component([])
    elif mode == 'face':
        ax.face.hide_component([])        
    elif mode == 'edge':
        ax.edge.hide_component([])                
    elif mode == 'point':
        ax.point.hide_component([])                        
    else:
        pass
    viewer.draw_all()    

def hide_elem(evt):
    viewer = evt.GetEventObject().GetTopLevelParent()
    mode = viewer._sel_mode

    ax = viewer.get_axes()
    if mode == 'volume':
        facesa = []
        facesb = []        
        s, v = viewer._s_v_loop['geom']
        selected_volume = viewer._dom_bdr_sel                
        for key in v.keys():
            if key in selected_volume:
                facesa.extend(v[key])
            else:
                facesb.extend(v[key])
        facesa = np.unique(np.array(facesa))
        facesb = np.unique(np.array(facesb))
        new_hide = list(np.setdiff1d(facesa, facesb, True))
        idx = ax.face.hidden_component
        idx = list(set(idx+new_hide))
        ax.face.hide_component(idx)        
    elif mode == 'face':
        idx = ax.face.getSelectedIndex()
        idx = list(set(ax.face.hidden_component+idx))        
        ax.face.hide_component(idx)        
    elif mode == 'edge':
        idx = ax.edge.getSelectedIndex()
        idx = list(set(ax.edge.hidden_component+idx))        
        ax.edge.hide_component(idx)                
    elif mode == 'point':
        pass
    else:
        pass
    viewer.canvas.unselect_all()
    viewer.draw_all()
    
btask = [('gdot',    fdot,  2, 'select dot', select_dot),
         ('gedge',   fedge, 2, 'select edge', select_edge),
         ('gface',   fface, 2, 'select face', select_face),
         ('gdomain', fdom,  2, 'select domain', select_volume),
         ('---', None, None, None),
         ('gshow',   show,  0, 'show all', show_all),
         ('ghide',   hide,  0, 'hide selection', hide_elem),]         
            
