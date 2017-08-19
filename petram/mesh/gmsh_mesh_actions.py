import numpy as np

import petram.debug as debug
dprint1, dprint2, dprint3 = debug.init_dprints('GmshMeshActions')

from petram.phys.vtable import VtableElement, Vtable
from petram.mesh.gmsh_mesh_model import GmshMeshActionBase



class MeshData(object):
    def __init__(self, lines, num_entities):
        self.lines = lines
        
        self.done = {"point":[],     #0D element
                     "edge": [],     #1D element
                     "face": [],     #2D element
                     "volume": []}   #3D element
        self.num_entities = {"point": num_entities[0],
                             "edge": num_entities[1],
                             "face": num_entities[2],
                             "volume": num_entities[3],}
        
        
    def append(self, c):
        self.lines.append(c)

    def show_hide_gid(self, gid, mode = "edge"):
        if gid.strip() == "":
            return "" # dont do anything
        elif gid == "*":
            show(self, gid, mode = mode)            
        elif gid == "remaining":
            if self.done[mode] == "*": return "" # already all done
            show(self, "*", mode = mode)
            hide(self, ','.join([str(x) for x in self.done[mode]]),
                 mode = mode)            
        else:
            hide(self, "*", mode = mode)
            show(self, gid, mode = mode)
            
        if gid == "*":        
            self.done[mode] = "*"
        elif gid == "remaining":
            gid = self.get_remaining_txt(mode)
            self.done[mode] = "*"            
        else:
            gidnum = [int(x) for x in gid.split(',')]
            for x in gidnum:
                if not x in self.done[mode]:
                    self.done[mode].append(x)
        return gid
                    
    def get_remaining_txt(self, mode = "edge"):
        if self.done[mode] == "*": return ''
        ll = [x+1 for x in range(self.num_entities[mode])]
        for x in self.done[mode]:ll.remove(x)
        if len(ll) == 0:
            self.done[mode] = "*"
            return ''
        else:
            return ','.join([str(x) for x in ll])
            
        

data = (('geom_id', VtableElement('geom_id', type='string',
                                   guilabel = 'Line#',
                                   default = "remaining", 
                                   tip = "Line ID" )),
        ('num_seg', VtableElement('radius', type='float',
                                   guilabel = 'Number of segments',
                                   default = 5, 
                                   tip = "Number of segments" )),
        ('progression', VtableElement('progression', type='float',
                                   guilabel = 'Progression',
                                   default = 0, 
                                   tip = "Progression" )),
        ('bump', VtableElement('bump', type='float',
                                   guilabel = 'Bump',
                                   default = 0, 
                                   tip = "Bump" )),)

        
    
class TransfiniteLine(GmshMeshActionBase):
    vt = Vtable(data)    
    def build_mesh(self, lines):
        gid, nseg, p, b = self.vt.make_value_or_expression(self)
        gid = lines.show_hide_gid(gid, mode="edge")
        if gid == "": return
        transfinite(lines, gid, mode = 'Line', nseg=nseg,
                    progression = p,  bump = b)
        mesh(lines, dim = 1)            


    def get_element_selection(self):
        self.vt.preprocess_params(self)                
        ret, mode = super(TransfiniteLine, self).get_element_selection()
        try:
            ret['edge'] = [int(x) for x in self.geom_id.split(',')]
        except:
            pass
        return ret, 'edge'
