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
        ax.faces.hide_component([])
    elif mode == 'face':
        ax.faces.hide_component([])        
    elif mode == 'edge':
        ax.edges.hide_component([])                
    elif mode == 'point':
        ax.points.hide_component([])                        
    else:
        pass
    viewer.draw_all()    

def hide_elem(evt):
    viewer = evt.GetEventObject().GetTopLevelParent()
    mode = viewer._sel_mode

    ax = viewer.get_axes()
    if mode == 'volume':
        idx = ax.faces.getSelectedIndex()
        idx = list(set(ax.faces.hidden_component+idx))
        #### I should select faces which only belongs to
        #### a particular volume !!!
        ax.faces.hide_component(idx)        
    elif mode == 'face':
        idx = ax.faces.getSelectedIndex()
        idx = list(set(ax.faces.hidden_component+idx))        
        ax.faces.hide_component(idx)        
    elif mode == 'edge':
        idx = ax.edges.getSelectedIndex()
        idx = list(set(ax.edges.hidden_component+idx))        
        ax.edges.hide_component([])                
    elif mode == 'point':
        pass
    else:
        pass
    viewer.draw_all()
            
btask = [('dot',    fdot,  2, 'select dot', select_dot),
         ('edge',   fedge, 2, 'select edge', select_edge),
         ('face',   fface, 2, 'select face', select_face),
         ('domain', fdom,  2, 'select domain', select_volume),
         ('---', None, None, None),
         ('show',   show,  0, 'show all', show_all),
         ('hide',   hide,  0, 'hide selection', hide_elem),]         
            
