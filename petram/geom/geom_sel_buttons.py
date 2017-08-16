from petram.utils import get_pkg_datafile
import petram.geom

fdot = get_pkg_datafile(petram.pi, 'icon',  'dot.png')
fedge = get_pkg_datafile(petram.pi, 'icon', 'line.png')
fface = get_pkg_datafile(petram.pi, 'icon', 'face.png')
fdom = get_pkg_datafile(petram.pi, 'icon', 'domain.png')

def select_dot(evt):
    print(evt.GetEventObject())
    
def select_edge(evt):
    print(evt.GetEventObject())
    
def select_face(evt):
    print(evt.GetEventObject())
    
def select_dom(evt):
    print(evt.GetEventObject())    
    
btask = [('dot',    fdot,  2, 'select dot', select_dot),
         ('edge',   fedge, 2, 'select edge', select_edge),
         ('face',   fface, 2, 'select face', select_face),
         ('domain', fdom,  2, 'select domain', select_dom),]
            
