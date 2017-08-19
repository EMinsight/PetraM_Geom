import sys

def show(gid, mode = "Line"):
    lines= []
    if gid == "*":
        lines.append('Show "*";')
    elif gid == "": return
    else:
        gid = [str(x) for x in gid.split(',')]
        txt = 'Recursive Show {{ {}; }}'.format(
                mode + '{{{}}}'.format(','.join(gid)))
        lines.append(txt)
    return lines

def hide(gid, mode = "Line"):
    lines = []    
    if gid == "*":
        lines.append('Hide "*";')
    elif gid == "": return        
    else:
        gid = [str(x) for x in gid.split(',')]
        if len(gid) == 0: return
        txt = 'Recursive Hide {{ {}; }}'.format(
                mode + '{{{}}}'.format(','.join(gid)))
        lines.append(txt)
    return lines
        
def mesh(dim = 1):
    lines = []
    lines.append('Mesh.MeshOnlyVisible=1;')
    lines.append('Mesh ' + str(dim) + ';')
    return lines

def transfinite(gid, mode = 'Line', nseg='',
                progression = 0, bump = 0):
    lines = []
    c = 'Transfinite '+mode
    if gid == "*":
        c += ' "*" = '
        print(c)
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

class Mesher(object):
    def __init__(self, CharacteristicLengthMax = 1e20,
                       CharacteristicLengthMin = 1):


        self.clmax = CharacteristicLengthMax
        self.clmin = CharacteristicLengthMin

        self.sequence = []

    def add_xxx(self, name, *gids, **kwargs):
        m = getattr(self, 'add_' + name)
        self.sequence.append([name, gids, kwargs])

    def add_transfinite(self, *gids, **kwargs):
        self.sequence.append(['transfinite', gids, kwargs])

    def generate(self):
        lines = []
        thismodule = sys.modules[__name__]
        for proc, gids, kwargs in self.sequence:
            f = getattr(thismodule, proc+'1d')
            lines.append(f(*gids, **kwargs))
            f = getattr(thismodule, proc+'2d')
            lines.append(f(*gids, **kwargs))
            f = getattr(thismodule, proc+'3d')
            lines.append(f(*gids, **kwargs))
        return lines
    
            
            
