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
        ptx, l, s = geom_data

    l1 = np.unique(np.hstack([s[k] for k in src]).flatten())
    l2 = np.unique(np.hstack([s[k] for k in dst]).flatten())
    p1p = np.unique(np.hstack([l[k] for k in l1]).flatten())
    p2p = np.unique(np.hstack([l[k] for k in l2]).flatten())

    p1 = ptx[p1p-1]
    p2 = ptx[p2p-1]
    n1 = normal2points(p1)
    n2 = normal2points(p2)

    ax = np.cross(n1, n2)
    an = np.arcsin(np.linalg.norm(ax))    
    ax = ax/np.linalg.norm(ax)
    
    def find_mapping(p1, p3):
        # try all transpose
        for i in range(len(p1)):
            d = p3[0]- p1[i]
            p3t = p3 - d
            mind = np.array([np.min(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1])
            if np.all(mind < mind_eps): 
                 mapping = [np.argmin(np.sqrt(np.sum((p3t - pp)**2,1))) for pp in p1]
                 trans = d
                 return d, mapping
        return None, None
    
    if an*180./np.pi > min_angle:
        R = rotation_mat(ax, -an)
    else:
        R = np.diag([1,1,1.])
        an = 0.0

    # check two possible orientation        
    p3 = np.dot(R, p2.transpose()).transpose()
    d, mapping = find_mapping(p1, p3)
    if d is None:
        if an*180./np.pi > min_angle:
            R = rotation_mat(ax, np.pi-an)
        else:
            R = np.diag([1,1,1.])
            an = 0.0
        p3 = np.dot(R, p2.transpose()).transpose()
        d, mapping = find_mapping(p1, p3)
        if d is None:        
            assert False, "auto trans failed (no mapping between vertices)"
    p_pairs = dict(zip(p1p, p2p[mapping]))  #point mapping

    l2dict = {tuple(sorted(l[ll])):ll for ll in l2}
    l_pairs = {ll:l2dict[tuple(sorted((p_pairs[l[ll][0]],p_pairs[l[ll][1]])))] for ll in l1}

    affine = np.zeros((4,4), dtype=float)
    affine[:3,:3] = np.linalg.inv(R)
    affine[:3,-1] = np.dot(np.linalg.inv(R), d)
    affine[-1,-1] = 1.0

    px = np.dot(np.linalg.pinv(-R+np.diag((1,1,1))),-d)
    print("px, d", px, d)
    return ax, an, px, d, affine, p_pairs, l_pairs
