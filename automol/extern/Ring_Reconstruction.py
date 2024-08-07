""" 
Reconstruct ring geometry from puckering parameters
Code from Chan et al. 2021 https://doi.org/10.1021/acs.jcim.0c01144
"""
##################################
#### Ring Reconstruction (RR) ####
##################################


import numpy
from ..geom import (
                        translate_to_ring_center,
                        normal_to_ring_plane
                        )
# import py_rdl
# import Ring_Analysis as RA
# import bondtable

###############################
########## Fix Zero ###########
###############################

def fixzero(x):
    """
    Fixes precision issues
    """
    x_ = numpy.array([0.0]) if numpy.allclose(0,x, rtol=1e-06, atol=1e-08) else x
    return x_.item()


###################################################
##### Ring Reconstruction for MONOCYCLIC ring #####
###################################################
# def GetCoordinate(mol, ringpath):
#     """
#     Get molecule coordinates

#     Input:

#     mol: rdmol

#     ringpath: list

#     Return:

#     coordinate: list
#     """
#     coordinate = [list(mol.GetConformer().GetAtomPosition(x)) for x in ringpath]
#     return coordinate

# def GetRingBondLength(mol, ringpath):
#     """
#     Get bond length of the ring bonds

#     Input:

#     mol: rdmol

#     ringidx: list

#     Return:

#     bondlength: list

#     """
#     N = len(ringpath)
#     ringbond = [[ringpath[i], ringpath[(i+1)%N]] for i in range(N)]
#     molconf = mol.GetConformer()
#     bondlength = [rdMolTransforms.GetBondLength(molconf, *b) for b in ringbond]
#     return bondlength

# def SetInitialRingBondLength(mol, ringpath):
#     """
#     Set inital ring bond length for the reference planar ring

#     Input:

#     mol: rdmol

#     ringidx: list

#     Return:

#     bondlength: list

#     """
#     N = len(ringpath)
#     pair = [[ringpath[i], ringpath[(i+1)%N]] for i in range(N)]
#     first_ele = [mol.GetAtomWithIdx(x[0]).GetAtomicNum() for x in pair]
#     second_ele = [mol.GetAtomWithIdx(x[1]).GetAtomicNum() for x in pair]
#     sorted_ele = [sorted(p) for p in zip(first_ele,second_ele)]
#     bond = [int(mol.GetBondBetweenAtoms(*x).GetBondType()) for x in pair]
#     ele_bond = [(sorted_ele[i][0],sorted_ele[i][1], bond[i]
#     ) for i,j in enumerate(bond)]
#     bondlength = [bondtable.BONDLENGTH_REF.get(eb,
#     "Unknown bond length, please update BOND LENGTH table"
#     ) for eb in ele_bond]
#     return bondlength

# def GetRingBondAng(mol, ringpath):
#     """
#     Get bond angles of the ring

#     Input:

#     mol: rdmol

#     ringpath: list

#     Return:

#     bondang: list (output in radian)

#     """
#     N = len(ringpath)
#     atoms = [[ringpath[i], ringpath[(i+1)%N], ringpath[(i+2)%N]] for i in range(N)]
#     molconf = mol.GetConformer()
#     bondang =[rdMolTransforms.GetAngleRad(molconf, x[0], x[1], x[2]) for x in atoms]
#     return bondang

# def SetInitialRingBondAng(mol, ringpath):
#     """
#     Get bond angles of the ring

#     Input:

#     mol: rdmol

#     ringpath: list

#     Return:

#     bondang: list (output in radian)

#     """
#     N = len(ringpath)
#     atoms = [[ringpath[i], ringpath[(i+1)%N], ringpath[(i+2)%N]] for i in range(N)]
#     first_ele = [mol.GetAtomWithIdx(x[0]).GetAtomicNum() for x in atoms]
#     second_ele = [mol.GetAtomWithIdx(x[1]).GetAtomicNum() for x in atoms]
#     third_ele = [mol.GetAtomWithIdx(x[2]).GetAtomicNum() for x in atoms]
#     sorted_ele = list(zip(first_ele, second_ele, third_ele))
#     bond = []
#     for a in atoms:
#         bond.append([int(mol.GetBondBetweenAtoms(a[i],a[i+1]
#                                                  ).GetBondType()) for i in range(2)])
#     elements_and_bond_order = [sorted_ele[i]+tuple(bond[i]
#                               ) for i,j in enumerate(bond)]
#     bondang = [bondtable.BONDANGLE_REF.get(x,
#             "Unknown bond length, please update BOND ANGLE table"
#             ) for x in elements_and_bond_order]
#     return bondang

def GetZ(ringsize, q, ang):
    """
    Get Z position from desired puckering amplitude (q) and angle (ang). 
    Angle is in radians

    :param ringsize: ring size
    :type ringsize: int
    :param q: q pucker params
    :type q: list of floats
    :param ang: phi pucker params
    :type ang: list of floats
    :rtype Z: list of floats
    """
    qsize = len(q)
    angsize = len(ang)
    assert ((qsize+angsize)==(ringsize-3)) and 4<ringsize<16
    Z = []
    K = int((ringsize-1)/2)+1 if ringsize%2==1 else int(ringsize/2)
    if ringsize%2==1: # odd ring size
        for j in range(ringsize):
            tmp = numpy.sum([q[m-2]*numpy.cos(ang[m-2]
                  +2*numpy.pi*m*j/ringsize) for m in range(2,K)])
            Z.append(numpy.sqrt(2/ringsize)*tmp)
    else:  # even ring size
        for j in range(ringsize):
            tmp = numpy.sqrt(2/ringsize)*numpy.sum([q[m-2]*numpy.cos(
                ang[m-2]+2*numpy.pi*m*j/ringsize) for m in range(2,K
                )])+numpy.sqrt(1/ringsize)*q[-1]*(-1)**j
            Z.append(tmp)
    return Z

def UpdateBondLength(ringsize, Z, init_bondl):
    """
    Update bond length given new z position 

    :param ringsize: ring size
    :type ringsize: int
    :param Z: displacement of ring atoms from mean ring plane
    :type Z: list of floats
    :param init_bondl: ring bond lengths
    :type init_bondl: list of floats
    :rtype nr: list of floats
    """
    assert len(Z)==ringsize and len(init_bondl)==ringsize, "Inappropriate Input"
    N = ringsize
    nr = [numpy.sqrt(numpy.square(init_bondl[i])-numpy.square(
        Z[(i+1)%N]-Z[i])) for i in range(N)] #Order of bonds from 0-1 to last-0
    return nr

def UpdateBeta(ringsize, Z, r, beta):
    """
    Update bond angle given new z position

    :param ringsize: ring size
    :type ringsize: int
    :param Z: displacement of ring atoms from mean ring plane
    :type Z: list of floats
    :param r: ring bond lengths
    :type r: list of floats
    :param beta: ring angle amplitudes
    :type beta: list of floats
    """
    assert len(Z)==ringsize and len(r)==ringsize and len(beta)==ringsize
    N = ringsize
    r_new = UpdateBondLength(ringsize,Z,r)
    idxlist = [[i,(i+1)%N,(i+2)%N] for i in range(ringsize)]
    beta_new = []
    for i in idxlist:
        tmp_a = numpy.square(Z[i[2]]-Z[i[0]])-numpy.square(
            Z[i[1]]-Z[i[0]])-numpy.square(Z[i[2]]-Z[i[1]])
        tmp_b = 2*r[i[0]]*r[i[1]]*numpy.cos(beta[i[0]])
        tmp_c = 2*r_new[i[0]]*r_new[i[1]]
        if (tmp_a+tmp_b)/tmp_c > 0:
            beta_new.append(numpy.arccos(min(1.,(tmp_a+tmp_b)/tmp_c)))
        else:
            beta_new.append(numpy.arccos(max(-1.,(tmp_a+tmp_b)/tmp_c)))
    return beta_new

def RotationMatrix(phi):
    """
    Rotation matrix
    
    :param phi: rotation angle
    :type phi: float
    :rtype rotation: numpy.array
    """
    rotation = numpy.array([(numpy.cos(phi), -numpy.sin(phi)
                             ), (numpy.sin(phi),numpy.cos(phi))])
    return rotation

def SegCoord(segment, r, beta):
    """
    Coordinates of each segment.

    :param segment: ring atoms indexes
    :type segment: list of ints
    :param r: ring bond lengths
    :type r: list of floats
    :param beta: ring angle amplitudes
    :type beta: list of floats
    :rtype coordinate: list
    """
    coordinate = []
    segsize = len(segment)
    segment_r = [r[x] for x in segment]
    segment_beta = [beta[x] for index, x in enumerate(segment)]
    alpha = []
    gamma = []
    coordinate = [numpy.array((0,0)),numpy.array((segment_r[0],0))]
    slength = [segment_r[0]]
    if segsize>=3:
        for index in range(segsize-2):
            if index==0:
                x = slength[-1] + r[index+1]*numpy.cos(numpy.pi-beta[index])
                y = r[index+1]*numpy.sin(numpy.pi-beta[index])
                coordinate.append(numpy.array((x,y)))
            else:
                x = slength[-1] + r[index+1]*numpy.cos(numpy.pi-alpha[index-1])
                y = r[index+1]*numpy.sin(numpy.pi-alpha[index-1])
                rota = RotationMatrix(gamma[index-1])
                coord = numpy.matmul(rota,numpy.array((x,y)))
                coordinate.append(coord)
            slength.append(numpy.linalg.norm(coordinate[-1]))
            alpha.append(segment_beta[index+1]-numpy.arcsin(numpy.sin(
                segment_beta[index])*slength[index]/slength[index+1]))
            gamma.append(numpy.arctan2(y,x))
    return coordinate

def RingPartition(ringsize, z, r, beta):
    """
    Initialize coordinates for three segments.

    :param ringsize: ring size
    :type ringsize: int
    :param z: displacement of ring atoms from mean ring plane
    :type z: list of floats
    :param r: ring bond lengths
    :type r: list of floats
    :param beta: ring angle amplitudes
    :type beta: list of floats
    :rtype finalcoord: numpy.array
    """
    # divide the ring into segments and initialise the coordinate of the segments
    location = list(range(ringsize))
    assert 5<=ringsize<=16
    if 5<=ringsize<=7:
        segment1 = location[0:3]
        segment2 = location[2:-1]
        segment2.reverse()
        segment3 = location[-2:]+[location[0]]
        segment3.reverse()
    elif 8<=ringsize<=10:
        segment1 = location[0:4]
        segment2 = location[3:-2]
        segment2.reverse()
        segment3 = location[-3:]+[location[0]]
        segment3.reverse()
    elif 11<=ringsize<=13:
        segment1 = location[0:5]
        segment2 = location[4:-3]
        segment2.reverse()
        segment3 = location[-4:]+[location[0]]
        segment3.reverse()
    else: #ringsize<=16 and ringsize>=14:
        segment1 = location[0:6]
        segment2 = location[5:-4]
        segment2.reverse()
        segment3 = location[-5:] + [location[0]]
        segment3.reverse()
    segcoord_1_init = SegCoord(segment1, r, beta)
    segcoord_2_init = SegCoord(segment2, r, beta)
    segcoord_3_init = SegCoord(segment3, r, beta)
    Reflection = numpy.array((-1,1))
    OPsq = numpy.inner(segcoord_1_init[-1], segcoord_1_init[-1])
    PQsq = numpy.inner(segcoord_2_init[-1], segcoord_2_init[-1])
    OQsq = numpy.inner(segcoord_3_init[-1], segcoord_3_init[-1])
    segcoord_1 = [Reflection*item for item in segcoord_1_init]
    segcoord_2 = [x + numpy.sqrt((OQsq,0)) for x in segcoord_2_init]
    segcoord_3 = [numpy.array(x) for x in segcoord_3_init]
    # Link segment together
    xp = (OPsq+OQsq-PQsq)/(2*numpy.sqrt(OQsq))
    yp = numpy.sqrt(OPsq-numpy.square(xp))
    phi1 = numpy.arctan2(segcoord_1[-1][1],segcoord_1[-1][0])
    phi2 = numpy.arctan2(segcoord_2[-1][1], segcoord_2[-1][0]-numpy.sqrt(OQsq))
    phi3 = numpy.arctan2(segcoord_3[-1][1], segcoord_3[-1][0])
    phiseg1 = numpy.arctan2(yp,xp)
    phiseg2 = numpy.arctan2(yp,xp-numpy.sqrt(OQsq))
    sigma1 = numpy.abs(phi1-phiseg1)
    sigma2 = numpy.abs(phiseg2-phi2)
    sigma3 = numpy.abs(phi3)
    Rsigma1 = RotationMatrix(-sigma1)
    Rsigma2 = RotationMatrix(sigma2)
    Rsigma3 = RotationMatrix(-sigma3)
    coordinate_1 = [numpy.array((0,0))]
    seg1_size = len(segcoord_1)
    for i in range(1,seg1_size-1):
        coordinate_1.append(numpy.matmul(Rsigma1,segcoord_1[i]))
    coordinate_1.append(numpy.array((xp,yp)))
    #### Check Here ####
    coordinate_2 = []
    seg2_size = len(segcoord_2)
    for i in range(seg2_size-2,0,-1):
        tmp = numpy.sqrt((OQsq,0))
        coordinate_2.append(tmp + numpy.matmul(Rsigma2, (segcoord_2[i]-tmp)))
    coordinate_3 = [numpy.sqrt((OQsq,0))]
    seg3_size = len(segcoord_3)
    for i in range(seg3_size-2,0,-1):
        coordinate_3.append(numpy.matmul(Rsigma3, segcoord_3[i]))
    coordinate = coordinate_1 + coordinate_2 + coordinate_3
    Rg = numpy.sum(coordinate,axis=0)
    phig = numpy.arctan2(Rg[1],Rg[0]) + numpy.pi/2
    Rphig = RotationMatrix(-phig)
    newcoord = [numpy.matmul(Rphig, coordinate[i]-Rg).tolist(
                )+[z[i]] for i in range(ringsize)]
    origin = numpy.mean(newcoord,axis=0)
    finalcoord = numpy.array(newcoord)-origin
    return finalcoord

####################################################
#### Call this function to reconstruct the ring ####
####################################################
#def SetRingPuckerCoords(mol, ringpath, amplitude, angle, init_bondl, init_bondang):
def SetRingPuckerCoords(ringpath, amplitude, angle, init_bondl, init_bondang):
    """
    Reconstruct the ring from given puckering amplitude and puckering angle

    :param ringpath: ring atoms
    :type ringpath: list of ints
    :param amplitude: q pucker params
    :type amplitude: list of floats
    :param angle: phi pucker params
    :type angle: list of floats
    :param init_bondl: bond lengths between ring atoms
    :type init_bondl: list of floats
    :param init_bondang: angles between ring atoms
    :type init_bondang: list of floats
    :rtype newcoord: list of lists
    :rtype newcoord: list of floats    
    """
#    molcenter = numpy.array(GetCoordinate(mol, ringpath)).mean(axis=0)
    N = len(ringpath)
    newZ = GetZ(N, amplitude, angle)
    new_bondl = UpdateBondLength(N, newZ, init_bondl)
    new_bondang = UpdateBeta(N, newZ, init_bondl, init_bondang)
    newcoord = numpy.array(RingPartition(N, newZ, new_bondl, new_bondang))
    numpy.set_printoptions(formatter={'float_kind':'{:f}'.format})
    return newcoord.tolist(),newZ



########################################
########## Ring Substituents  ##########
########################################

def GetRingSubstituentPosition(coord, sub_coord,ring_sub_idx):
    """
    Determine alpha and beta paramters for a given substituent on a ring.

    :param coord: ring cartesian coordinates
    :type coord: list of tuples (x,y,z)
    :param sub_coord: coordinates of substituent
    :type sub_coord: tuple (x,y,z)
    :param ring_sub_idx: index of ring atom attached to sub
    :type ring_sub_idx: int 
    :rtype alpha: float (0, numpy.pi)
    :rtype beta: float  (-numpy.pi,*numpy.pi )
    """
    ring_coord = numpy.array(coord)
    substituent_coord = numpy.array([coord[ring_sub_idx],sub_coord])
    alpha=0
    beta=0
    ring_coord_ = translate_to_ring_center(ring_coord)
    substituent_coord_ = substituent_coord - ring_coord.mean(axis=0)
    S = numpy.diff(substituent_coord_,axis=0)
    s = S/numpy.linalg.norm(S)
    n = normal_to_ring_plane(ring_coord_)
    alpha = numpy.arccos(fixzero(numpy.dot(s,n)))
    alpha = alpha.item()
    R = numpy.array(substituent_coord_)[0]
    U = R - numpy.dot(R,n)*n
    u = U/numpy.linalg.norm(U)
    v = numpy.cross(n,u)
    su = fixzero(numpy.dot(s,u))
    sv = fixzero(numpy.dot(s,v))
    beta = numpy.arctan2(-sv,su)
    beta = beta.item()
    return alpha, beta

def SetRingSubstituentPosition(coord, alpha, beta, sub_coord, bondlength):
    """
    Update ring subtituent position. Bond length is fixed.

    :param coord: ring cartesian coordinates
    :type coord: list of tuples (x,y,z)
    :param alpha: first parameter for ring subs
    :type alpha: float (0, numpy.pi)
    :param beta: second parameter for ring subs
    :type beta: float  (-numpy.pi,*numpy.pi )
    :param sub_coord: coordinates of ring atom to which sub is attached
    :type sub_coord: tuple (x,y,z)
    :param bondlength: length of bond (ring_atom - sub)
    :type bondlength: float 
    :rtype ring_substituent_pos: tuple (x,y,z)
    """
    ring_coord = numpy.array(coord)
    ring_at_coord = numpy.array(sub_coord)
    ring_coord_ = translate_to_ring_center(ring_coord)
    n = normal_to_ring_plane(ring_coord_)
    # ring atom position WITH RESPECT TO RING CENTER
    R = numpy.array(ring_at_coord) - ring_coord.mean(axis=0)
    U = R - numpy.dot(R,n)*n # u is projection of R orthogonal to n
    u = U/numpy.linalg.norm(U)
    v = numpy.cross(n,u) # v is orthogonal to n and u
    x = fixzero(bondlength*numpy.sin(alpha)*numpy.cos(-beta))
    y = fixzero(bondlength*numpy.sin(alpha)*numpy.sin(-beta))
    z = fixzero(bondlength*numpy.cos(alpha))
    b = numpy.array([x,y,z])
    T = numpy.array([u,v,n]).T
    ring_substituent_pos = numpy.matmul(T,b) + ring_at_coord
    return tuple(ring_substituent_pos.tolist())
