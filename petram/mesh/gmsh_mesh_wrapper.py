from __future__ import print_function
import numpy as np
import gmsh

gmsh_init = False

debug = True
def dprint(*args):
    if debug:
        print(*args)
    
class GMSHMeshWrapper(object):
    workspace_base = 1
    def __init__(self, format=2.2):
        global gmsh_init
        if not gmsh_init:
           gmsh.initialize()
        gmsh.clear()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", format)
        gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)
        gmsh_init = True
        self.add_model('main1')
        
    def load_brep(self, filename):
        self.switch_model('main1')
        self.target = filename
        
        gmsh.model.occ.importShapes(filename)
        gmsh.model.occ.synchronize()

        self.geom_info = self.read_geom_info()
        return self.current
    
    def read_geom_info(self):
        from petram.geom.read_gmsh import read_loops2
        self.hide_all()        
        gmsh.model.mesh.generate(1)

        #ptx, p, l, s, v = read_loops2(gmsh)
        return read_loops2(gmsh)
        
    def add_model(self, name):
        gmsh.model.add(name)
        name = gmsh.model.list()[-1]
        self.current = name
        return name
        
    def switch_model(self, name = 'main1'):
        gmsh.model.setCurrent(name)
        self.current = name
        return name
    
    def prep_workspace(self):
        name = 'ws'+str(GMSHMeshWrapper.workspace_base)
        GMSHMeshWrapper.workspace_base += 1
        
        self.add_model(name)
        gmsh.model.setCurrent(name)        
        gmsh.model.occ.importShapes(self.target)
        self.hide_all()        
        gmsh.model.mesh.generate(1)
        gmsh.model.occ.synchronize()
        
        return name
        
    def show_only(self, dimtags, recursive=False):
        self.hide_all()
        gmsh.model.setVisibility(dimtags, True, recursive = recursive)
        
    def hide(self, dimtags, recursive=False):
        gmsh.model.setVisibility(dimtags, False, recursive = recursive)
        
    def hide_all(self):
        ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, False)
        
    def show_all(self):
        ent = gmsh.model.getEntities()
        gmsh.model.setVisibility(ent, True)

    def delete_all_except(self, dim, tags):
        for d in [3, 2, 1]:
            ent = gmsh.model.getEntities(d)            
            if dim == d:
                ent = [(d, t) for d, t in ent if not t in tags]
                gmsh.model.occ.remove(ent, recursive=True)
                break
            else:
                gmsh.model.occ.remove(ent)
        gmsh.model.occ.synchronize()        

    def transfer_mesh(self, dimtags, ws1, ws2, resursive = True):
        '''
        transfer mesh data of dimtags from ws1 to ws2
        '''
        if resursive:
            dimtags = self.expand_dimtags(dimtags)

        old_model = self.current
        
        data = []
        self.switch_model(ws1)
        for dim, tag in dimtags:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            data.append((dim, tag, ndata, edata))
        self.switch_model(ws2)            
        for dim, tag, ndata, edata in data:
            gmsh.model.mesh.setNodes(dim, tag, ndata[0], ndata[1], ndata[2])
            gmsh.model.mesh.setElements(dim, tag, edata[0], edata[1], edata[2])

        self.switch_model(old_model)
            
    def save_mesh(self, filename, ws=''):
        '''
        save mesh
        '''
        current = self.current
        self.switch_model(ws)
        gmsh.write(filename)
        self.switch_model(self.current)
        
    def save_mesh_all(self, filename_base):
        '''
        save mesh for all models
        '''
        workspaces = gmsh.model.list()
        for ws in workspaces:
            f = filename_base+'_'+ws+'.msh'
            self.save_mesh(f, ws=ws)

    def expand_dimtags(self, dimtags):
        '''
        expand all subdims
        '''
        volumes = []
        surfaces = []
        lines = []
        points = []
        
        for dim, tag in dimtags:
            if dim == 3: volumes.append(tag)
            if dim == 2: surfaces.append(tag)
            if dim == 1: lines.append(tag)
            if dim == 0: points.append(tag)            
        for v in volumes:           
            surfaces.extend(self.geom_info[3][v])

        surfaces = list(set(surfaces))
        for s in surfaces:
            lines.extend(self.geom_info[2][s])
                        
        lines = list(set(lines))
        for l in lines:
            points.extend(self.geom_info[1][l])

        points = list(set(points))

        ans = ( [(0, p) for p in points] +
                [(1, l) for l in lines] +
                [(2, s) for s in surfaces] +
                [(3, v) for v in volumes] )
        return ans

    def list_entities(self):
        for x in gmsh.model.list():
            gmsh.model.setCurrent(x)
            print(x, gmsh.model.getEntities())
        gmsh.model.setCurrent(self.current)

    '''
    High-level interface test codes
    '''
    def extrude_surface_0D(self, vtag, tag1, tag2, nlayers):
        pass
    def extrude_surface_1D(self, vtag, tag1, tag2, nlayers):
        pass
    def extrude_surface_2D(self, vtag, tag1, tag2, nlayers):
        pass
        
    '''
    High-level interface test codes
    '''
    def extrude_surface_test(self, vtag, tag1, tag2, nlayers):
        
        from petram.geom.geom_utils import find_translate_between_surface
        from petram.geom.geom_utils import map_points_in_geom_info
        from petram.geom.geom_utils import map_lines_in_geom_info
        from petram.geom.geom_utils import map_surfaces_in_geom_info
        from petram.geom.geom_utils import map_volumes_in_geom_info        
        
        ptx, p, l, s, v = self.geom_info
        geom_data = (ptx, l, s, None)
        ax, an, px, d, affine, p_pairs, l_pairs = find_translate_between_surface(
                                                            [tag1], [tag2],
                                                            geom_data=geom_data)
        
        ents = gmsh.model.getBoundary(((2,tag1),), oriented=False)
        ents = list(ents)+[(2,tag1),]
        mdata = []
        for dim, tag in ents:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))
        
        ents = gmsh.model.getBoundary(((2,tag2),), oriented=False)        
        
        ws = self.prep_workspace()

        self.delete_all_except(2, [tag1])
        ret = gmsh.model.occ.extrude(((2,tag1),), d[0], d[1], d[2],
                                     numElements=[nlayers])
        gmsh.model.occ.synchronize()        

        info1 = self.geom_info
        info2 = self.read_geom_info()

        pmap, pmap_r = map_points_in_geom_info(info1, info2)
        lmap, lmap_r = map_lines_in_geom_info(info1, info2, pmap_r)
        smap, smap_r = map_surfaces_in_geom_info(info1, info2, lmap_r)
        vmap, vmap_r = map_volumes_in_geom_info(info1, info2, smap_r)
        dprint("point map", pmap)
        dprint("line map", lmap)

        self.show_all()
        gmsh.model.mesh.generate(1)
        ents = gmsh.model.getBoundary(((2,tag1),ret[0]), oriented=False)
        for dim, tag in ents:
            gmsh.model.mesh.setNodes(dim, tag, [], [], [])
            gmsh.model.mesh.setElements(dim, tag, [], [], [])

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        #print(node_map)
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # copy 1D elements on the source surface
        for d in mdata:
            dim, tag, ndata, edata = d 
            if dim == 2: break
            tag = lmap[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
        print("node map here", node_map1)
        # copy 1D elements on the destination surface
        l_pairs2 = {x:lmap[l_pairs[x]] for x in l_pairs}
        R = affine[:3,:3]
        D = affine[:3,-1]
        node_map2 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        for p in p_pairs:
            node_map2[info1[1][p]+1] = info2[1][pmap[p_pairs[p]]]+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim == 2: break
            tag = l_pairs2[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map2[i] = j

            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]
            #print("setting", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        # copy source mesh
        for d in mdata:
            dim, tag, ndata, edata = d            
            if dim == 1: continue
            
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j

            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)
            
        # copy dest mesh
        for d in mdata:
            dim, tag, ndata, edata = d            
            if dim == 1: continue
            tag = ret[0][1]
            
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map2[i] = j

            pos = np.array(pos).reshape(-1,3).transpose()
            pos = (np.dot(R, pos).transpose() + D).flatten()
            gmsh.model.mesh.setNodes(dim, tag, ntag2, pos, ndata[2])

            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map2[x] for x in item] for item in nodes]
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.hide(((2,tag1),ret[0]),)
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.generate(3)                    


        #
        #  copy back mash info
        #
        # first gather info
        ents = gmsh.model.getEntities()
        ents2 = gmsh.model.getBoundary(((2,tag1),), oriented=False)
        ents = [x for x in ents if not x in ents2 and x[0] != 0]
        ents = [x for x in ents if x[0] != 2 or x[1] != tag1]
        mdata = []
        for dim, tag in ents:
            ndata = gmsh.model.mesh.getNodes(dim, tag)
            edata = gmsh.model.mesh.getElements(dim, tag)
            mdata.append((dim, tag, ndata, edata))


        
        # collect data for node mapping
        ents2_1D =  gmsh.model.getEntities(1)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents2_1D]
        ntags_ws_1D = sum([x[0] for x in tmp], [])
        pos_ws_1D = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        ents2_2D = gmsh.model.getEntities(2)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents2_2D]
        ntags_ws_2D = sum([x[0] for x in tmp], [])
        pos_ws_2D = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)

        #node_map1 = {info2[1][k]+1:info1[1][pmap_r[k]]+1 for k in pmap_r}
        
        self.switch_model('main1')
        ents2 = [x[1] for x in gmsh.model.getBoundary(((3,vtag),))]
        ents2 = sum([[y[1] for y in
                      gmsh.model.getBoundary(((2,x),), oriented=False)]
                     for x in ents2],[])
        ents2 = list(set(ents2))
        ents_1D = [(1, x) for x in ents2]
        
        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}

        '''
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]
        print("nodemap1", node_map1)
        '''        
        #gmsh.model.occ.synchronize()
        #self.list_entities()                
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1

        # send back 1D data
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 1: continue
            print("copy from ", dim, tag)                        
            tag = lmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            #print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}
        
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]

        # send back 2D data
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 2: continue
            print("copy from ", dim, tag)            
            tag = smap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.show_all()            
        self.hide(((3, 2),), recursive=True)
        gmsh.model.mesh.generate(2)

        node_map1 = {info1[1][k]+1:info2[1][pmap[k]]+1 for k in pmap}
        node_map1 = {node_map1[x]:x for x in node_map1}
        
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents_1D]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_1D]
        for i, nt in zip(idx, ntags_ws_1D):
            node_map1[nt] = ntags_m[i]

        ents2 = gmsh.model.getBoundary(((3,vtag),), oriented=False)
        tmp = [gmsh.model.mesh.getNodes(dim, tag) for dim, tag in ents2]
        ntags_m = sum([x[0] for x in tmp], [])
        pos_m = np.array(sum([x[1] for x in tmp], [])).reshape(-1,3)
        idx = [np.argmin(np.sum((pos_m-p)**2, 1)) for p in pos_ws_2D]
        for i, nt in zip(idx, ntags_ws_2D):
            node_map1[nt] = ntags_m[i]
        noffset = max(gmsh.model.mesh.getNodes()[0])+1
        eoffset = max(sum(gmsh.model.mesh.getElements()[1],[]))+1
        
        for d in mdata:
            dim, tag, ndata, edata = d
            if dim != 3: continue
            print("copy from ", dim, tag)            
            tag = vmap_r[tag]
            ntag, pos, ppos = ndata
            ntag2 = range(noffset, noffset+len(ntag))
            noffset = noffset+len(ntag)
            for i, j in zip(ntag, ntag2): node_map1[i] = j
            #print("setting nodes", dim, tag, ntag2, ndata[1], ndata[2])
            gmsh.model.mesh.setNodes(dim, tag, ntag2, ndata[1], ndata[2])
            etypes, etags, nodes = edata

            etags2 = [range(eoffset, eoffset+len(etags[0]))]
            eoffset = eoffset+len(etags[0])
            nodes2 = [[node_map1[x] for x in item] for item in nodes]
            #print("setting elements", dim, tag, etypes, etags2, nodes2)
            gmsh.model.mesh.setElements(dim, tag, etypes, etags2, nodes2)

        self.show_all()            
        self.hide(((3, 2),), recursive=True)
        gmsh.model.mesh.generate(3)
            
            
