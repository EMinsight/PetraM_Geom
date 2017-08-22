import sys

def show(gid, mode = "Line", recursive = False):
    lines= []
    if gid == "*":
        lines.append('Show "*";')
    elif gid == "": return []
    else:
        gid = [str(x) for x in gid.split(',')]
        txt = 'Show {{ {}; }}'.format(
                mode + '{{{}}}'.format(','.join(gid)))
        if recursive: txt = 'Recursive ' + txt
        lines.append(txt)
    return lines

def hide(gid, mode = "Line", recursive = False):
    lines = []    
    if gid == "*":
        lines.append('Hide "*";')
    elif gid == "": return []    
    else:
        gid = [str(x) for x in gid.split(',')]
        if len(gid) == 0: return
        txt = 'Hide {{ {}; }}'.format(
                mode + '{{{}}}'.format(','.join(gid)))
        if recursive: txt = 'Recursive ' + txt        
        lines.append(txt)
    return lines
        
def mesh(dim = 1):
    lines = []
    lines.append('Mesh.MeshOnlyVisible=1;')
    lines.append('Mesh ' + str(dim) + ';')
    return lines

def transfinite(gid, mode = 'Line', nseg='',
                progression = 0, bump = 0, meshdim = 1):
    lines = []
    c = 'Transfinite '+mode
    if gid == "*":
        c += ' "*" = '
    else:
        gid = [str(x) for x in gid.split(',')]
        c += '{{ {} }}'.format(','.join(gid)) + ' = '
    c += str(nseg)
    if bump != 0:
        c += " Using Bump " + str(bump)
    if progression != 0:
        c += " Using Progression " + str(bump)    
    lines.append(c+';')
    return lines

def freemesh(gid, clmax = None, clmin = None):
    lines = []
    if clmax > 0:
        lines.append('Mesh.CharacteristicLengthMax = '+str(clmax) + ';')
    if clmin > 0:
        lines.append('Mesh.CharacteristicLengthMin = '+str(clmin) + ';')
    if len(lines) > 0:
        lines.append('Mesh.CharacteristicLengthExtendFromBoundary = 0;')
    return lines

def characteristiclength(gid, cl = 1e20):
    c = "Characteristic Length "
    gid = [str(x) for x in gid.split(',')]    
    c += '{{ {} }}'.format(','.join(gid)) + ' = ' +  str(cl) + ";"
    return [c,]

def embed(gid, embed_s=None, embed_l=None, embed_p=None):
    if ((embed_s is None) and (embed_l is None) and
        (embed_p is None)): return []
   
    gid = [str(x) for x in gid.split(',')]
    if len(gid) != 1:
        assert False, "embed destination should be one element"

    lines = []
    if embed_s is not None:
        sid = [str(x) for x in embed_s.split(',') if x]
        if sid:
            lines.append('Surface {{ {} }}'.format(','.join(sid)) +
                      ' In Volume {{ {} }}'.format(','.join(gid))+ ';')
    if embed_l is not None:
        lid = [str(x) for x in embed_l.split(',') if x]
        print(lid)
        if lid:        
            lines.append('Line {{ {} }}'.format(','.join(lid)) +
                      ' In Surface {{ {} }}'.format(','.join(gid))+ ';')
    if embed_p is not None:
        pid = [str(x) for x in embed_p.split(',') if x]
        if pid:        
            lines.append('Point {{ {} }}'.format(','.join(pid)) +
                      ' In Surface {{ {} }}'.format(','.join(gid))+ ';')
    return lines           

class GmshMesher(object):
    def __init__(self, num_entities,
                       CharacteristicLengthMax = 1e20,
                       CharacteristicLengthMin = 1):


        self.clmax = CharacteristicLengthMax
        self.clmin = CharacteristicLengthMin

        self.sequence = []
        self.done = {"Vertex":[],     #0D element
                     "Line": [],     #1D element
                     "Surface": [],     #2D element
                     "Volume": []}   #3D element
        self.num_entities = {"Vertex": num_entities[0],
                             "Line": num_entities[1],
                             "Surface": num_entities[2],
                             "Volume": num_entities[3],}
        
    def transfinite(self, gid, mode = 'Line', nseg='',
                    progression = 0, bump = 0, meshdim = 1, **kwargs):
        lines = []
        if meshdim == 0:
            if gid == 'remaining':
                gid = self.get_remaining_txt(mode = 'Line')
            lines.extend(transfinite(gid, mode = mode,
                         nseg = nseg,
                         progression = progression,
                                    bump = bump))
        else:
            if mode == 'Line' and meshdim != 1: return []
            if mode == 'Surface' and meshdim != 2: return []
        
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
        self.record_finished(gid, mode = mode)
        return lines

    def freemesh(self, gid,  mode='Line', clmax = -1, clmin = -1,
                 meshdim = 1, embed_s = None, embed_l = None, embed_p=None):
        '''
        freemesh  = unstructured volume/surface/line
        '''
        clmax  = self.clmax if clmax == -1 else clmax
        clmin  = self.clmin if clmin == -1 else clmin
        
        lines = []
        if meshdim == 3 and mode == 'Volume':
            pass
        elif meshdim == 2 and mode == 'Surface':
            pass
        elif meshdim == 1 and mode == 'edge':
            pass
        else:
            lines = embed(gid, embed_s=embed_s, embed_l=embed_l,
                          embed_p=embed_p)
            return lines

        x = self.show_hide_gid(gid, mode = mode)
        if len(x) == 0: return lines
        lines.extend(x)
        lines.extend(freemesh(gid, clmax=clmax, clmin=clmin))
        self.record_finished(gid, mode = mode)        
        return lines
    
    def characteristiclength(self, gid, cl = 1.0, meshdim = 0):
        if meshdim == 0:
            return characteristiclength(gid, cl = cl)
        else:
            return []
        
    def add(self, name, *gids, **kwargs):
        self.sequence.append([name, gids, kwargs])
        
    def reset_done(self):
        self.done = {"Vertex":[],     #0D element
                     "Line": [],     #1D element
                     "Surface": [],     #2D element
                     "Volume": []}   #3D element

    def generate(self):
        lines = []
        thismodule = sys.modules[__name__]
        for mdim in [0, 1, 2, 3]:
            for proc, gids, kwargs in self.sequence:
                f = getattr(self, proc)
                kwargs['meshdim'] = mdim
                x = f(*gids, **kwargs)
                if len(x) > 0:
                    lines.extend(x)
                    if mdim > 0:
                        lines.extend(mesh(dim = mdim))
            self.reset_done()                    
        lines.extend(hide("*", mode = "Volume"))
        return lines

    def show_hide_gid(self, gid, mode = "Line", recursive = True):
        lines = []
        if gid.strip() == "":
            return lines # dont do anything
        elif gid == "*":
            lines.extend(show(gid, mode = mode, recursive = recursive))
        elif gid == "remaining":
            if self.done[mode] == "*": return lines # already all done
            lines.extend(show("*", mode = mode, recursive=recursive))
            lines.extend(hide(','.join([str(x) for x in self.done[mode]]),
                              mode = mode, recursive=recursive))
        else:
            lines.extend(hide("*", mode = mode, recursive=recursive))
            lines.extend(show(gid, mode = mode, recursive=recursive))
        return lines

    def record_finished(self, gid, mode = 'Line'):
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
                    
    def get_remaining_txt(self, mode = "Line"):
        if self.done[mode] == "*": return ''
        ll = [x+1 for x in range(self.num_entities[mode])]
        for x in self.done[mode]:ll.remove(x)
        if len(ll) == 0:
            self.done[mode] = "*"
            return ''
        else:
            return ','.join([str(x) for x in ll])
    
            
            