import numpy as np
import scipy

def rotation_mat(ax, an):
    '''
    matrix to rotate arounc ax by an [rad]
    '''
    c = np.cos(an); s = np.sin(an)
    ax = ax/np.sqrt(np.sum(ax**2))
    R = np.array(
            [[c + (1-c)*ax[0]**2, ax[0]*ax[1]*(1-c)-ax[2]*s, ax[0]*ax[2]*(1-c)+ax[1]*s],
             [ax[0]*ax[1]*(1-c)+ax[2]*s, c + (1-c)*ax[1]**2,  ax[1]*ax[2]*(1-c)-ax[0]*s],
             [ax[0]*ax[2]*(1-c)-ax[1]*s, ax[1]*ax[2]*(1-c)+ax[0]*s, c + (1-c)*ax[2]**2]]
            )

    return R

def normal2points(p1, eps = 1e-13):
    '''
    normal vector defined by a surface made from the group of
    points
  
    p1 : [#, 3] matrix, where # is the number of points
    '''
    p1 = np.hstack([p1, np.atleast_2d(np.array([1]*len(p1))).transpose()])
    u, s, vh = np.linalg.svd(p1)
 
    null_mask = (s <= eps)
    if sum(null_mask) == 0:
        print("no null space??", p1, s)
    null_space = scipy.compress(null_mask, vh, axis=0)
    norm = null_space[0, :3]
    norm = norm / np.sqrt(np.sum(norm**2))
    return norm

def map_points_in_geom_info(info1, info2, th = 1e-15):
    '''
    info = ptx, l, s, v
        pts = array(:, 3)
        p   = point -> point index
        l   = line -> point
        s   = surface -> line
        v   = volume -> surface
    '''
    ptx1 = info1[0]
    ptx2 = info2[0]
    
    dist = np.array([np.min(np.sum((ptx1 - p)**2, 1))for p in ptx2])
    if np.any(dist > th):
        assert False, "could not able to find vertex mapping"


    #iverts -> p
    iv2p1 = {info1[1][k]:k    for k in info1[1]}
    iv2p2 = {info2[1][k]:k    for k in info2[1]}
    
    # {point in info2 : point in info1}
    pmap_r = {iv2p2[k]: iv2p1[np.argmin(np.sum((ptx1 - p)**2, 1))]
              for k,  p in enumerate(ptx2)}
    # {point in info1 : point in info2}    
    pmap = {pmap_r[k]:k   for k in pmap_r}

    return pmap, pmap_r

def map_lines_in_geom_info(info1, info2, pmap_r):
    lmap = {}
    lmap_r = {}

    for l in info2[2]:
        p1, p2 = pmap_r[info2[2][l][0]], pmap_r[info2[2][l][1]]
        for x in info1[2]:
            if (info1[2][x][0] == p1 and
                info1[2][x][1] == p2):
                lmap[x] = l
                lmap_r[l] = x
                break
            elif (info1[2][x][0] == p2 and
                  info1[2][x][1] == p1):
                lmap[x] = -l
                lmap_r[l] = -x
                break
            else:
                pass
        else:
            assert False, "could not find line mapping for "+str(l)
    return lmap, lmap_r

def map_surfaces_in_geom_info(info1, info2, lmap_r):
    smap = {}
    smap_r = {}

    for s in info2[3]:
        tmp = sorted([abs(lmap_r[x]) for x in info2[3][s]])
        for x in info1[3]:
            if sorted(info1[3][x]) == tmp:
                smap[x] = s
                smap_r[s] = x
                break
            else:
                pass
        else:
            assert False, "could not find surface mapping for "+str(s)
    return smap, smap_r

def map_volumes_in_geom_info(info1, info2, smap_r):
    vmap = {}
    vmap_r = {}

    for v in info2[4]:
        tmp = sorted([smap_r[x] for x in info2[4][v]])
        for x in info1[4]:
            if sorted(info1[4][x]) == tmp:
                vmap[x] = v
                vmap_r[v] = x
                break
            else:
                pass
        else:
            assert False, "could not find volume mapping for "+str(s)
    return vmap, vmap_r

                
def find_translate_between_surface(src, dst, geom=None,
                                   geom_data = None,
                                   min_angle = 0.1,
                                   mind_eps = 1e-10,
                                   axan = None):
    
    if geom is not None:
        ptx, cells, cell_data, l, s, v, geom = geom._gmsh4_data
    else:
        ptx, l, s, cell_data = geom_data

    l1 = np.unique(np.hstack([s[k] for k in src]).flatten())
    l2 = np.unique(np.hstack([s[k] for k in dst]).flatten())
    p1p = np.unique(np.hstack([l[k] for k in l1]).flatten())
    p2p = np.unique(np.hstack([l[k] for k in l2]).flatten())

    if cell_data is None:
        i1 = p1p-1
        i2 = p2p-1
    else:
        i1 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p1p]).flatten()
        i2 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p2p]).flatten()
    p1 = ptx[i1,:]
    p2 = ptx[i2,:]
    n1 = normal2points(p1)
    n2 = normal2points(p2)

    if axan is None:
        ax = n1
        an = 0.0


        '''
        ax = ax/np.linalg.norm(ax)    
        n3 = np.cross(ax, n1)
        xx = np.sum(n2*n1)
        yy = np.sum(n2*n3)
        #an = np.arcsin(np.linalg.norm(ax))
        an = np.arctan2(yy, xx)
        #print("p2, axis angle", xx, yy, p2, ax, an)        
        '''
    else:
        ax, an = axan
        ax = np.array(ax, dtype=float)
        ax = ax/np.linalg.norm(ax)            
        an = np.pi/180.*an
        

    def find_mapping(ax, an, p1, p2):
        if an != 0.0:
            R = rotation_mat(ax, -an)
        else:
            R = np.diag([1,1,1.])

    # check two possible orientation        
        p3 = np.dot(R, p2.transpose()).transpose()
        
        # try all transpose
        #print("p1, p3 (1)", p1, p3)        
        for i in range(len(p1)):
            d = p3[0]- p1[i]
            p3t = p3 - d
            mind = np.array([np.min(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1])
            if np.all(mind < mind_eps): 
                 mapping = [np.argmin(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1]
                 trans = d
                 return d, mapping, R
        return None, None, None
    
    if abs(an*180./np.pi) < min_angle:
        an = 0.0
    if abs(abs(an*180./np.pi)-180.) < min_angle:
        an = 0.0
        
    d, mapping, R = find_mapping(ax, an, p1, p2)

    if d is None:
        if an > 0.:
            an = -np.pi + an
        else:
            an = np.pi + an
        d, mapping, R = find_mapping(ax, an, p1, p2)        
        if d is None:        
            assert False, "auto trans failed (no mapping between vertices)"

    p_pairs = dict(zip(p1p, p2p[mapping]))  #point mapping

    #print("l1", [(ll, l[ll]) for ll in l1])
    #print("l2", [(ll, l[ll]) for ll in l2])    
    l2dict = {tuple(sorted(l[ll])):ll for ll in l2}

    l_pairs = {ll:l2dict[tuple(sorted((p_pairs[l[ll][0]],p_pairs[l[ll][1]])))] for ll in l1}

    affine = np.zeros((4,4), dtype=float)
    affine[:3,:3] = np.linalg.inv(R)
    affine[:3,-1] = np.dot(np.linalg.inv(R), d)
    affine[-1,-1] = 1.0

    px = np.dot(np.linalg.pinv(-R+np.diag((1,1,1))),-d)
    #print("px, d", px, d)
    return ax, an, px, d, affine, p_pairs, l_pairs

def find_rotation_between_surface(src, dst, geom=None,
                                   geom_data = None,
                                   min_angle = 0.1,
                                   mind_eps = 1e-10,
                                   axan = None):
    
    if geom is not None:
        ptx, cells, cell_data, l, s, v, geom = geom._gmsh4_data
    else:
        ptx, l, s, cell_data = geom_data

    l1 = np.unique(np.hstack([s[k] for k in src]).flatten())
    l2 = np.unique(np.hstack([s[k] for k in dst]).flatten())
    p1p = np.unique(np.hstack([l[k] for k in l1]).flatten())
    p2p = np.unique(np.hstack([l[k] for k in l2]).flatten())

    if cell_data is None:
        i1 = p1p-1
        i2 = p2p-1
    else:
        i1 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p1p]).flatten()
        i2 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p2p]).flatten()
    p1 = ptx[i1,:]
    p2 = ptx[i2,:]
    n1 = normal2points(p1)
    n2 = normal2points(p2)

    #print(p1, p2, n1, n2, axan)
    if axan is None:
        M = np.vstack((n1, n2))
        b = np.array([np.sum(n1*p1[0]), np.sum(n2*p2[0])])
        
        from scipy.linalg import null_space
        from numpy.linalg import lstsq


        ax =null_space(M).flatten()
        px, res, rank, s = lstsq(M, b, rcond=None)
        
        pp1 = px - p1[0] - np.sum((px-p1[0])*ax)*ax
        pp2 = px - p2[0] - np.sum((px-p2[0])*ax)*ax

        pp1 = pp1/np.linalg.norm(pp1)
        pp2 = pp2/np.linalg.norm(pp2)

        s = np.mean(np.cross(pp1, pp2)[ax!=0]/ax[ax!=0])
        c = np.sum(pp1 * pp2)
        #xx = np.sum(n2*n1)
        #yy = np.sum(n2*n3)
        #an = np.arcsin(np.linalg.norm(ax))
        an = np.arctan2(s, c)

        
        #print("p2, axis angle", px, ax, an)        
    else:
        ax, an = axan
        ax = np.array(ax, dtype=float)
        ax = ax/np.linalg.norm(ax)            
        an = np.pi/180.*an
        

    def find_mapping(ax, an, p1, p2):
        if an != 0.0:
            R = rotation_mat(ax, -an)
        else:
            R = np.diag([1,1,1.])

    # check two possible orientation        
        p3 = np.dot(R, p2.transpose()).transpose()
        
        # try all transpose
        #print("p1, p3 (1)", p1, p3)        
        for i in range(len(p1)):
            d = p3[0]- p1[i]
            p3t = p3 - d
            mind = np.array([np.min(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1])
            if np.all(mind < mind_eps): 
                 mapping = [np.argmin(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1]
                 trans = d
                 return d, mapping, R
        return None, None, None
    
    if abs(an*180./np.pi) < min_angle:
        an = 0.0
    if abs(abs(an*180./np.pi)-180.) < min_angle:
        an = 0.0
        
    d, mapping, R = find_mapping(ax, an, p1, p2)

    if d is None:
        if an > 0.:
            an = -np.pi + an
        else:
            an = np.pi + an
        d, mapping, R = find_mapping(ax, an, p1, p2)        
        if d is None:        
            assert False, "auto trans failed (no mapping between vertices)"

    p_pairs = dict(zip(p1p, p2p[mapping]))  #point mapping

    #print("l1", [(ll, l[ll]) for ll in l1])
    #print("l2", [(ll, l[ll]) for ll in l2])    
    l2dict = {tuple(sorted(l[ll])):ll for ll in l2}

    l_pairs = {ll:l2dict[tuple(sorted((p_pairs[l[ll][0]],p_pairs[l[ll][1]])))] for ll in l1}


    affine = np.zeros((4,4), dtype=float)
    affine[:3,:3] = np.linalg.inv(R)
    affine[:3,-1] = np.dot(np.linalg.inv(R), d)
    affine[-1,-1] = 1.0
    
    if axan is None:
        px = np.dot(np.linalg.pinv(-R+np.diag((1,1,1))),-d)
    #print("px, d", px, d)
    return ax, an, px, d, affine, p_pairs, l_pairs
