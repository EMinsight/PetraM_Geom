try:
   import gmsh
   has_gmsh = True   
except ImportError:
   has_gmsh = False
   assert False, "gmsh api is not found"

class GmshMesher(object):
    def __init__(self, geom_root,
                       CharacteristicLengthMax = 1e20,
                       CharacteristicLengthMin = 1,
                       MeshAlgorithm = "Automatic",
                       MeshAlgorithm3D = "Delaunay"):

        cells = geom_root._gmsh4_data[1]        
        self.geom = geom_root._gmsh4_data[-1]
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
        
        self.num_entities = {"Vertex": len(geom.model.getEntities(0)),
                             "Line":   len(geom.model.getEntities(1)),
                             "Surface":len(geom.model.getEntities(2)),
                             "Volume": len(geom.model.getEntities(3)) }

        l = geom_root._gmsh4_data[3]
        s = geom_root._gmsh4_data[4]
        v = geom_root._gmsh4_data[5]                
        self.entity_relations = {"Volume": v, 
                                 "Surface": s, 
                                 "Line":  l }

        self.ietg = 0

    def add(self, name, *gids, **kwargs):
        self.sequence.append([name, gids, kwargs])

    def new_etg(self):
        self.ietg = self.ietg + 1
        return 'etg'+str(self.ietg)
    
    def reset_done(self):
        self.done = {"Vertex":[],     #0D element
                     "Line": [],     #1D element
                     "Surface": [],     #2D element
                     "Volume": []}   #3D element

    def generate(self):
        lines = []
        thismodule = sys.modules[__name__]
        gmsh.option.setNumber("Mesh.Algorithm", self.algorithm)
        gmsh.option.setNumber("Mesh.Algorithm3D", self.algorithm3d)
        
        if self.clmax > 0:
            gmsh.option.setNumber("Mesh.CharacteristicLengthMax",
                                  self.clmax)
        if self.clmin > 0:
            gmsh.option.setNumber("Mesh.CharacteristicLengthMin",
                                  self.clmin)
        
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
        return lines, max_mdim
        
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
