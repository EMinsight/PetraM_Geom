import numpy as np
import scipy

def rotation_mat(ax, an):
    '''
    matrix to rotate arounc ax by an [rad]
    '''
    c = np.cos(an); s = np.sin(an)
    R = np.array(
            [[c + (1-c)*ax[0]**2, ax[0]*ax[1]*(1-c)-ax[2]*s, ax[0]*ax[2]*(1-c)+ax[1]*s],
             [ax[0]*ax[1]*(1-c)+ax[2]*s, c + (1-c)*ax[1]**2,  ax[1]*ax[2]*(1-c)-ax[0]*s],
             [ax[0]*ax[2]*(1-c)-ax[1]*s, ax[1]*ax[2]*(1-c)+ax[0]*s, c + (1-c)*ax[2]**2]]
            )

    return R

def normal2points(p1, eps = 1e-15):
    '''
    normal vector defined by a surface made from the group of
    points
  
    p1 : [#, 3] matrix, where # is the number of points
    '''
    p1 = np.hstack([p1, np.atleast_2d(np.array([1]*len(p1))).transpose()])
    u, s, vh = np.linalg.svd(p1)

    null_mask = (s <= eps)
    if sum(null_mask) == 0:
        print("no null space??", p1)
    null_space = scipy.compress(null_mask, vh, axis=0)
    norm = null_space[0, :3]
    norm = norm / np.sqrt(np.sum(norm**2))
    return norm

def find_translate_between_surface(src, dst, geom=None,
                   geom_data = None,
                   min_angle = 0.1, 
                   mind_eps = 1e-10):
    
    if geom is not None:
        ptx, cells, cell_data, l, s, v, geom = geom._gmsh4_data
    else:
        ptx, l, s, cell_data = geom_data

    l1 = np.unique(np.hstack([s[k] for k in src]).flatten())
    l2 = np.unique(np.hstack([s[k] for k in dst]).flatten())
    p1p = np.unique(np.hstack([l[k] for k in l1]).flatten())
    p2p = np.unique(np.hstack([l[k] for k in l2]).flatten())

    i1 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p1p]).flatten()
    i2 = np.array([np.where(cell_data['vertex']['geometrical'] == ii)[0] for ii in p2p]).flatten()
    p1 = ptx[i1,:]
    p2 = ptx[i2,:]
    n1 = normal2points(p1)
    n2 = normal2points(p2)

    ax = np.cross(n1, n2)
    ax = ax/np.linalg.norm(ax)    
    n3 = np.cross(ax, n1)
    xx = np.sum(n2*n1)
    yy = np.sum(n2*n3)
    #an = np.arcsin(np.linalg.norm(ax))
    an = np.arctan2(yy, xx)

    print("p2, axis angle", xx, yy, p2, ax, an)
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

    print("p_pairs here", p_pairs)
    print("l1", [(ll, l[ll]) for ll in l1])
    print("l2", [(ll, l[ll]) for ll in l2])    
    l2dict = {tuple(sorted(l[ll])):ll for ll in l2}
    l_pairs = {ll:l2dict[tuple(sorted((p_pairs[l[ll][0]],p_pairs[l[ll][1]])))] for ll in l1}

    affine = np.zeros((4,4), dtype=float)
    affine[:3,:3] = np.linalg.inv(R)
    affine[:3,-1] = np.dot(np.linalg.inv(R), d)
    affine[-1,-1] = 1.0

    px = np.dot(np.linalg.pinv(-R+np.diag((1,1,1))),-d)
    #print("px, d", px, d)
    return ax, an, px, d, affine, p_pairs, l_pairs
