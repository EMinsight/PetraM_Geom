from __future__ import print_function
import sys
import numpy as np
from collections import OrderedDict

MeshAlgorithm= OrderedDict((("MeshAdap", 1), ("Automatic", 2), ("Delaunay", 3),
                ("Frontal", 6), ("BAMG", 7), ("DelQuad", 8),
                ("default", 2)))

MeshAlgorithm3D= OrderedDict((("Delaunay",1), ("New Delaunay",2),
                              ("Frontal", 4), 
                              ("Frontal Hex", 6), ("MMG3D", 7),
                              ("R-tree", 9), ("default", 1)))
## note : 
##   Frontal Delaunay (3D) : deplicated.
        
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

def transfiniteL(gid, nseg='', progression = 0, bump = 0, meshdim = 1):
    lines = []
    mode = 'Line'
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
        c += " Using Progression " + str(progression)    
    lines.append(c+';')
    return lines

def transfiniteS(gid, points = None):
    if points is None: points = []
    lines = []
    mode = 'Surface'
    c = 'Transfinite '+mode
    if gid == "*":
        c += ' "*" = '
    else:
        gid = [str(x) for x in gid.split(',')]
        c += '{{ {} }}'.format(','.join(gid))
    if len(points) != 0:
        c += " = {"+','.join([str(x) for x in points])+"}"
    lines.append(c+';')
    return lines

def freemesh(gid, clmax=None, clmin=None):
    lines = []
    if clmax > 0:
        lines.append('Mesh.CharacteristicLengthMax = ' + str(clmax) + ';')
    if clmin > 0:
        lines.append('Mesh.CharacteristicLengthMin = ' + str(clmin) + ';')
    if len(lines) > 0:
        lines.append('Mesh.CharacteristicLengthExtendFromBoundary = 0;')
    else:
        lines.append('Mesh.CharacteristicLengthExtendFromBoundary = 1;')
    return lines

def characteristiclength(gid, cl = 1e20):
    c = "Characteristic Length "
    gid = [str(x) for x in gid.split(',')]    
    c += '{{ {} }}'.format(','.join(gid)) + ' = ' +  str(cl) + ";"
    return [c,]

def rotate(axis, x0, angle):
    taxis = [str(x) for x in axis]
    tx0 = [str(x) for x in x0]    
    return ('Rotate'+'{' + '{{ {} }}'.format(','.join(taxis)) + ','  
                         + '{{ {} }}'.format(','.join(tx0))  + ',' 
                         + str(angle) + '}')


def translate(delta):
    dd = [str(d) for d in delta]
    return 'Translate'+'{{ {} }}'.format(','.join(dd))    

def periodic(mode, gid, sid, transform):
    txt = 'Periodic '+ mode + ' '
    lines = []
    gid = [str(x) for x in gid.split(',')]
    sid = [str(x) for x in sid.split(',')]
    
    txt  += '{{ {} }}'.format(','.join(gid)) + ' = '
    txt  += '{{ {} }}'.format(','.join(sid))    
    txt  += ' ' + transform + ';'
    lines.append(txt)
    return lines

def boundary(mode, etg, gid):
    txt =  etg + "() = Unique(Abs(Boundary{" + mode
    txt +=  "{" + gid + "}"
    txt +=  ";}));"
    return [txt]

def embed(gid, embed_s="", embed_l="", embed_p=""):
    if ((embed_s is "") and (embed_l is "") and
        (embed_p is "")): return []

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
    def __init__(self, entities,
                       geom_coords,
                       CharacteristicLengthMax = 1e20,
                       CharacteristicLengthMin = 1,
                       MeshAlgorithm = "Automatic",
                       MeshAlgorithm3D = "Delaunay"):
        
        self.geom_coords = geom_coords
        self.clmax = CharacteristicLengthMax
        self.clmin = CharacteristicLengthMin
        self.algorithm = MeshAlgorithm
        self.algorithm3d = MeshAlgorithm3D

        self.sequence = []
        self.transform = {}
        self.done = {"Vertex":[],     #0D element
                     "Line": [],     #1D element
                     "Surface": [],     #2D element
                     "Volume": []}   #3D element
        num_entities = entities[0]
        self.num_entities = {"Vertex": num_entities[0],
                             "Line": num_entities[1],
                             "Surface": num_entities[2],
                             "Volume": num_entities[3],}
        self.entity_relations = entities[1]        
        self.ietg = 0
        
    def new_etg(self):
        self.ietg = self.ietg + 1
        return 'etg'+str(self.ietg)
    
    def transfinite_line(self, gid,  nseg='',
                    progression = 0, bump = 0, meshdim = 1, **kwargs):
        mode = "Line"
        lines = []
        if meshdim == 0:
            if gid == 'remaining':
                gid = self.get_remaining_txt(mode = 'Line')
            lines.extend(transfiniteL(gid,
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
    
    def transfinite_surface(self, gid, edges = None, meshdim=1, **kwargs):
        mode = 'Surface'
        lines = []
        if meshdim == 0:
            if gid == 'remaining':
                gid = self.get_remaining_txt(mode = 'Surface')
            # translate from edges to points
            rel = self.entity_relations
            if edges is not None:
                sid = int(gid)
                lids = rel['Surface'][sid]
                points = []
                print(edges)
                for e in edges:
                    if e.strip() == "": continue
                    eids = [int(x) for  x in e.split(',')]
                    print(e, eids)
                    for eid in eids:
                        assert (eid in lids), "edge is not part of surface: " + str(eid)
                    pts = list(np.hstack([rel['Line'][eid] for eid in eids]))
                    print(e, pts)                    
                    cc = []
                    for pt in pts:
                        if pts.count(pt) == 1: cc.append(pt)
                    if len(points) == 0:
                        points = cc
                    else:
                        if points[0] == cc[0]:
                            points = [cc[1]] + points
                        elif points[0] == cc[1]:
                            points = [cc[0]] + points
                        elif points[-1] == cc[0]:
                            points = points + [cc[1]]
                        elif points[-1] == cc[1]:
                            points = points + [cc[0]]
                        else:
                            pass
                print('points', points)                         
            else:
                points = None
            lines.extend(transfiniteS(gid, points = points))
        else:
            if mode == 'Line' and meshdim != 1: return []
            if mode == 'Surface' and meshdim != 2: return []
        
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
        self.record_finished(gid, mode = mode)
        return lines
    

    def freemesh(self, gid,  mode='Line', clmax=-1, clmin=-1,
                 meshdim=1, embed_s="", embed_l="", embed_p=""):
        '''
        freemesh  = unstructured volume/surface/line
        '''
        clmax  = self.clmax if clmax == -1 else clmax
        clmin  = self.clmin if clmin == -1 else clmin
        
        lines = []
        if meshdim == 3 and mode == 'Volume':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
        elif meshdim == 2 and mode == 'Volume':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
            if self.done['Surface']:
                lines.extend(hide(','.join([str(x) for x in self.done['Surface']]),
                                 mode = 'Surface'))
        elif meshdim == 1 and mode == 'Volume':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
            if self.done['Line']:
                lines.extend(hide(','.join([str(x) for x in self.done['Line']]),
                                 mode = 'Line'))
            if self.done['Surface']:
                lines.extend(hide(','.join([str(x) for x in self.done['Surface']]),
                                  mode = 'Surface', recursive = True))
        elif meshdim == 2 and mode == 'Surface':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
        elif meshdim == 1 and mode == 'Surface':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)
            if self.done['Line']:
                lines.extend(hide(','.join([str(x) for x in self.done['Line']]),
                                 mode = 'Line'))
        elif meshdim == 1 and mode == 'Line':
            x = self.show_hide_gid(gid, mode = mode)
            if len(x) == 0: return lines
            lines.extend(x)

        elif meshdim == 0:
            lines = embed(gid, embed_s=embed_s, embed_l=embed_l,
                          embed_p=embed_p)
            return lines
        else:
            return []

        lines.extend(freemesh(gid, clmax=clmax, clmin=clmin))
        self.record_finished(gid, mode = mode)                
        return lines

    def rotate(self, gid, src="", meshdim=0, transform=''):
        if meshdim!=0: return []
        gid = [int(x) for x in gid.split(',')]
        src = [int(x) for x in src.split(',')]
        if len(gid) != 1 or len(src) != 1:
            assert False, "source and dentination array length should be one"
        X, cells, pt_data, cell_data, field_data = self.geom_coords
        
        g = np.where(cell_data['triangle']['geometrical'] == gid[0])[0][0]
        s = np.where(cell_data['triangle']['geometrical'] == src[0])[0][0]

        g1 = X[cells['triangle'][g]]
        s1 = X[cells['triangle'][s]]

        norms = np.cross(s1[1]-s1[0], s1[2]-s1[0])  # normal on src
        normg = np.cross(g1[1]-g1[0], g1[2]-g1[0])  # normal on dst

        def norm(x): return np.sqrt(np.sum(x**2))

        norms = norms/norm(norms)
        normg = normg/norm(normg)

        #angle = np.arcsin(np.sum(norms*normg))
        axis = np.cross(norms, normg)
        axis = axis/norm(axis)
        
        m = np.vstack((norms, normg, axis))
        b = np.array([np.sum(norms*s1[0]), np.sum(normg*g1[0]), 0])
        x0 = np.dot(np.linalg.inv(m), b)

        d = sorted([(norm(np.cross(s1[k]-x0, axis)), s1[k] - x0)
                     for k in range(3)],
                     key=lambda xxxx: xxxx[0])
        p1 = d[-1][1]
        d = sorted([( norm(np.cross(g1[k]-x0, axis)), g1[k] - x0)
                     for k in range(3)],
                     key=lambda xxxx: xxxx[0])        
        p2 = d[-1][1]
        
        p1 = p1 - np.sum(p1*axis)*axis
        p2 = p2 - np.sum(p2*axis)*axis
        p1 = p1/norm(p1)
        p2 = p2/norm(p2)

        pos_p2 = np.cross(axis, p1)
        sign_angle = np.sum(pos_p2*p2)
        angle = np.arcsin(norm(np.cross(p1,p2)))
        if np.sum(p1*p2)<0: angle = np.pi - angle
        if sign_angle < 0: angle = - angle
        
        print('rotate', axis, x0, angle*180/np.pi)
        self.transform[transform] = ('rotate', axis, x0, angle,)
        return []
        
    def translate(self, gid, src="", meshdim=0, transform=''):
        if meshdim!=0: return []
        
        gid = [int(x) for x in gid.split(',')]
        src = [int(x) for x in src.split(',')]
        if len(gid) != 1 or len(src) != 1:
            assert False, "source and dentination array length should be one"
        X, cells, pt_data, cell_data, field_data = self.geom_coords
        
        g = np.where(cell_data['vertex']['geometrical'] == gid[0])[0]
        s = np.where(cell_data['vertex']['geometrical'] == src[0])[0]        
        a = X[cells['vertex'][g]]
        b = X[cells['vertex'][s]]
        self.transform[transform] = ('translate', (a-b).flatten())
        print('translate', (a - b).flatten())
        return []
        
    def copymesh(self, gid, src="", meshdim=0, transform='', mode='',
                       etg = None):
        etg1, etg2 = etg
        trans = self.transform[transform]
        f = getattr(sys.modules[__name__], trans[0])
        trans_txt = f(*trans[1:])
        
        x1 = self.show_hide_gid(gid, mode = mode)
        if not x1: return []
        x2 = self.show_hide_gid(src, mode = mode)
        if not x2: return []
        xx = self.show_hide_gid(','.join([gid, src]), mode = mode)        
        lines = []
        if meshdim == 0:
            lines.extend(boundary(mode, etg1, gid))
            lines.extend(boundary(mode, etg2, src))
            mode2 = 'Line' if mode == 'Surface' else 'Point'
            lines.extend(periodic(mode2, etg1+"()", etg2+"()", trans_txt))                
            lines.extend(periodic(mode,  gid, src, trans_txt))
        elif meshdim == 2 and mode == 'Surface':
            lines.extend(xx)
            self.record_finished(gid, mode = mode)
            self.record_finished(src, mode = mode)                        
        elif meshdim == 1 and mode == 'Surface':
            lines.extend(self.show_hide_gid(','.join([etg1+"()", etg2+"()"]),
                                            mode = 'Line'))
            self.record_finished(etg1+"()", mode = 'Line')
            self.record_finished(etg2+"()", mode = 'Line')
        elif meshdim == 1 and mode == 'Line':
            lines.extend(xx)
            self.record_finished(gid, mode = mode)
            self.record_finished(src, mode = mode)                        
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
        lines.append("Mesh.Algorithm = " +
                     str(MeshAlgorithm[self.algorithm])+';')
        lines.append("Mesh.Algorithm3D = " +
                     str(MeshAlgorithm3D[self.algorithm3d])+';')
        if self.clmax > 0:        
            lines.append('Mesh.CharacteristicLengthMax = ' +
                         str(self.clmax) + ';')
        if self.clmin > 0:                
            lines.append('Mesh.CharacteristicLengthMin = ' +
                         str(self.clmin) + ';')
        
        max_mdim = 0
        for mdim in [0, 1, 2, 3]:
            for proc, gids, kwargs in self.sequence:
                f = getattr(self, proc)
                kwargs['meshdim'] = mdim
                x = f(*gids, **kwargs)
                if len(x) > 0:
                    lines.extend(x)
                    if mdim > 0:
                        lines.extend(mesh(dim = mdim))
                        max_mdim = mdim
            self.reset_done()                    
        lines.extend(hide("*", mode = "Volume"))
        return lines, max_mdim

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
            try:
               gidnum = [int(x) for x in gid.split(',')]
            except:
               self.done[mode].append(gid)
               return
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
    
def write_entities_relations(geo_text, writev = -1, writes = -1):
    #from petram.geom.gmsh_geom_model import read_loops
    #geom_root = self.root()['Geometry'][self.geom_group]        
    #s, v = read_loops(geom_root._txt_unrolled)

    has_finalv = False
    has_finals = False
    has_finall = False
    for line in geo_text:
        if line.startswith('final_v[]'): has_finalv = True
        if line.startswith('final_s[]'): has_finals = True
        if line.startswith('final_l[]'): has_finall = True

    t1 = 'final_s() = Unique(Abs(Boundary{ Volume{final_v()}; }));'
    t2 = 'final_l() = Unique(Abs(Boundary{ Surface{final_s()}; }));'
    t3 = 'final_p() = Unique(Abs(Boundary{ Line{final_l()}; }));'

    tt1 = ['txt = "";',
           'For ii In {0 : #final_v[]-1}',
           '   pts() = Unique(Abs(Boundary{Volume{final_v[ii]}; }));',
           '   txt = StrCat("Boundary(Volume{", Sprintf("%g", final_v[ii]), "})=");',
           '   For iii In {0 : #pts[]-1}',
           '       txt=StrCat(txt, " ", Sprintf("%g", pts[iii]));',
           '   EndFor',
           '   Printf(txt);',
           'EndFor']
    tt2 = ['txt = "";',
           'For ii In {0 : #final_s[]-1}', 
           '   pts() = Unique(Abs(Boundary{Surface{final_s[ii]}; }));',
           '   txt = StrCat("Boundary(Surface{", Sprintf("%g", final_s[ii]), "})=");',
           '   For iii In {0 : #pts[]-1}',
           '       txt=StrCat(txt, " ", Sprintf("%g", pts[iii]));',
           '   EndFor',
           '   Printf(txt);',
           'EndFor']
    tt3 = ['txt = "";',
           'For ii In {0 : #final_l[]-1}',
           '   pts() = Unique(PointsOf{ Line{final_l[ii]}; });',
           '   txt = StrCat("PointsOf(Line{", Sprintf("%g", final_l[ii]), "})=");',
           '   For iii In {0 : #pts[]-1}',
           '       txt=StrCat(txt, " ", Sprintf("%g", pts[iii]));',
           '   EndFor',
           '   Printf(txt);',
           'EndFor',]


    ## if volume (surface loop) exitsts but there is only one volume,
    ## set it to final_v
    #len(v.keys()) == 1:
    #        ipv = str(v.keys()[0])
    if writev != -1:
        ipv = str(writev)
        geo_text.append('final_v[] = {'+ipv+'};')
        has_finalv = True
        has_finals = False
    ## if surface (line loop)  exitsts but there is only one volume,
    ## set it to final_s
    if writes != -1:    
    #if len(s.keys()) == 1:
        ips = str(writes)
        geo_text.append('final_s[] = {'+ips+'};')
        has_finals = True
    if has_finalv:
        geo_text.extend([t1, t2, t3]+ tt1 + tt2 + tt3)
    if has_finals:
        geo_text.extend([t2, t3] + tt2 + tt3)
    if has_finall:
        geo_text.extend([t3] + tt3)
    return geo_text

def write_physical(geo_text):
    has_finalv = False
    has_finals = False
    has_finall = False
    for line in geo_text:
        if line.startswith('final_v[]'):has_finalv = True
        if line.startswith('final_s[]'): has_finals = True
        if line.startswith('final_l[]'): has_finall = True
    if has_finalv:
        has_finals = False
        has_finall = False        
    if has_finals:
        has_finall = False        

    tt1= ['ipv = 0;',
          'For ii In {0 : #final_v[]-1}',
          '   ipv = ipv+1;',
          '   Physical Volume (StrCat("volume", Sprintf("%g", final_v[ii])),ipv) = {final_v[ii]};',
          'EndFor',]
    tt2 = ['ips = 0;',
           'For ii In {0 : #final_s[]-1}',
           '   ips = ips+1;',
           '   Physical Surface (StrCat("surface", Sprintf("%g", final_s[ii])),ips) = {final_s[ii]};',
           'EndFor']
    tt3 = ['ipl = 0;',
           'For ii In {0 : #final_l[]-1}',
           '   ipl = ipl+1;',
           '   Physical Line (StrCat("line", Sprintf("%g", final_l[ii])),ipl) = {final_l[ii]};',
          'EndFor',]

    if has_finalv:
        geo_text.extend(tt1 + tt2 + tt3)
    if has_finals:
        geo_text.extend(tt2 + tt3)
    if has_finall:
        geo_text.extend(tt3)
    return geo_text
