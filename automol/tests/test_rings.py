""" test ring functionality in automol.graph
"""
from automol import graph


def test__rings():
    """ test graph.rings
    """
    c5h5n5o_cgr = (
        {0: ('C', 1, None), 1: ('C', 0, None), 2: ('C', 0, None),
         3: ('C', 0, None), 4: ('C', 0, None), 5: ('N', 2, None),
         6: ('N', 0, None), 7: ('N', 0, None), 8: ('N', 0, None),
         9: ('N', 1, None), 10: ('O', 1, None)},
        {frozenset({10, 4}): (1, None), frozenset({8, 2}): (1, None),
         frozenset({0, 6}): (1, None), frozenset({9, 3}): (1, None),
         frozenset({1, 2}): (1, None), frozenset({3, 7}): (1, None),
         frozenset({2, 5}): (1, None), frozenset({1, 6}): (1, None),
         frozenset({0, 7}): (1, None), frozenset({9, 4}): (1, None),
         frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})

    assert graph.rings(c5h5n5o_cgr) == (
        ({0: ('C', 1, None), 1: ('C', 0, None), 3: ('C', 0, None),
          6: ('N', 0, None), 7: ('N', 0, None)},
         {frozenset({0, 6}): (1, None), frozenset({3, 7}): (1, None),
          frozenset({0, 7}): (1, None), frozenset({1, 6}): (1, None),
          frozenset({1, 3}): (1, None)}),
        ({1: ('C', 0, None), 2: ('C', 0, None), 3: ('C', 0, None),
          4: ('C', 0, None), 8: ('N', 0, None), 9: ('N', 1, None)},
         {frozenset({8, 2}): (1, None), frozenset({9, 3}): (1, None),
          frozenset({1, 2}): (1, None), frozenset({9, 4}): (1, None),
          frozenset({1, 3}): (1, None), frozenset({8, 4}): (1, None)})
    )


def test__ring_systems():
    """ test automol.graph.ring_systems
    """
    # molecule:
    # InChI=1S/C19H30/c1-2-4-14-10-12(9-13(14)3-1)5-7-17-16-8-6-15-11-
    # 18(16)19(15)17/h12-19H,1-11H2/
    gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 2, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
            6: ('C', 2, None), 7: ('C', 2, None), 8: ('C', 1, None),
            9: ('C', 2, None), 10: ('C', 2, None), 11: ('C', 2, None),
            12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 2, None),
            15: ('C', 2, None), 16: ('C', 2, None), 17: ('C', 2, None),
            18: ('C', 2, None)},
           {frozenset({9, 13}): (1, None), frozenset({3, 6}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({11, 12}): (1, None),
            frozenset({13, 14}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
            frozenset({12, 13}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({1, 7}): (1, None), frozenset({12, 15}): (1, None),
            frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
            frozenset({16, 15}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({16, 17}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({18, 4}): (1, None), frozenset({17, 14}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({18, 10}): (1, None),
            frozenset({8, 11}): (1, None)})
    rsys = graph.ring_systems(gra)
    assert len(rsys) == 2

    rsy_rngs = list(map(graph.rings, rsys))
    assert tuple(map(len, rsy_rngs)) == (3, 2)


def test__ring_systems_decomposed_atom_keys():
    """ test automol.graph.ring_systems_decomposed_atom_keys
    """
    # molecule:
    # InChI=1S/C19H30/c1-2-4-14-10-12(9-13(14)3-1)5-7-17-16-8-6-15-11-
    # 18(16)19(15)17/h12-19H,1-11H2/
    gra = ({0: ('C', 1, None), 1: ('C', 1, None), 2: ('C', 2, None),
            3: ('C', 1, None), 4: ('C', 1, None), 5: ('C', 1, None),
            6: ('C', 2, None), 7: ('C', 2, None), 8: ('C', 1, None),
            9: ('C', 2, None), 10: ('C', 2, None), 11: ('C', 2, None),
            12: ('C', 1, None), 13: ('C', 1, None), 14: ('C', 2, None),
            15: ('C', 2, None), 16: ('C', 2, None), 17: ('C', 2, None),
            18: ('C', 2, None)},
           {frozenset({9, 13}): (1, None), frozenset({3, 6}): (1, None),
            frozenset({0, 5}): (1, None), frozenset({11, 12}): (1, None),
            frozenset({13, 14}): (1, None), frozenset({3, 5}): (1, None),
            frozenset({0, 2}): (1, None), frozenset({1, 4}): (1, None),
            frozenset({12, 13}): (1, None), frozenset({0, 1}): (1, None),
            frozenset({1, 7}): (1, None), frozenset({12, 15}): (1, None),
            frozenset({6, 7}): (1, None), frozenset({8, 9}): (1, None),
            frozenset({16, 15}): (1, None), frozenset({4, 5}): (1, None),
            frozenset({16, 17}): (1, None), frozenset({2, 3}): (1, None),
            frozenset({18, 4}): (1, None), frozenset({17, 14}): (1, None),
            frozenset({8, 10}): (1, None), frozenset({18, 10}): (1, None),
            frozenset({8, 11}): (1, None)})

    decomps = graph.ring_systems_decomposed_atom_keys(gra)
    assert decomps == (((0, 1, 4, 5), (0, 2, 3, 5), (1, 7, 6, 3)),
                       ((8, 9, 13, 12, 11), (13, 14, 17, 16, 15, 12)))

# a1 = +/-q
# a2 = +/-a1


def test__ring_puckering():
    smi = 'CC1CCCCC1'
    import automol
    import numpy
    ich = automol.smiles.inchi(smi)
    # gra = automol.inchi.graph(ich)
    # zma = automol.graph.zmatrix(gra)
    geo = automol.inchi.geometry(ich)
    zma = automol.geom.zmatrix(geo)
    gra = automol.zmat.graph(zma)
    rings_atoms = automol.graph.rings_atom_keys(gra)
    val_dct = automol.zmat.value_dictionary(zma)
    coos = automol.zmat.coordinates(zma)
    geo = automol.zmat.geometry(zma)
    da_names = automol.zmat.dihedral_angle_names(zma)

    check_dct = {
        'dist': 3.5e-1,
        'coulomb': 1.5e-2,
       # 'stereo': None,
        'tors': None
    }
    for ring_atoms in rings_atoms:
        rotate_hyds = []
        ngbs = automol.graph.atom_sorted_neighbor_atom_keys(gra, ring_atoms[0])
        symbs = automol.geom.symbols(geo)
        for ngb in ngbs:
            if symbs[ngb] == 'H':
                rotate_hyds.append(ngb)
        ring_value_dct = {}
        for da in da_names:
            da_idxs = list(coos[da])[0]
            if len(list(set(da_idxs) & set(ring_atoms))) == 4:
                print(da, da_idxs)
                ring_value_dct[da] = val_dct[da]
        dist_value_dct = {}
        for i in range(len(ring_atoms)):
            dist_value_dct[i] = automol.zmat.distance(zma, ring_atoms[i-1], ring_atoms[i])

        samp_range_dct = {}
        for key, value in ring_value_dct.items():
           samp_range_dct[key] = (value - numpy.pi/4, value + numpy.pi/4)

        zmas = automol.zmat.samples(zma, 300, samp_range_dct)
        unique_geos = []
        unique_zmas = []
        for nzma in zmas:
            condition1 = True
            condition2 = False
            for i in range(len(ring_atoms)):
                if abs(dist_value_dct[i] - automol.zmat.distance(nzma, ring_atoms[i-1], ring_atoms[i])) > .3:
                    condition1 = False
            if condition1:
                condition2 = True
                new_geo = automol.zmat.geometry(nzma)
                for i in range(len(ring_atoms)):
                    angle_atoms = [ring_atoms[i], ring_atoms[i-1], ring_atoms[i-2]]
                    if automol.geom.central_angle(new_geo, *angle_atoms, degree=True) < 94.:
                        condition2 = False
            if condition2:    
                #if rotate_hyds:
                #    print('new\n')
                #    print(automol.geom.string(new_geo) + '\n')
                #    xyz_coords = automol.geom.coordinates(new_geo)
                #    rotate_angle = automol.geom.dihedral_angle(new_geo, *ring_atoms[-1:] + ring_atoms[:3])
                #    rotate_xyzs = tuple(a - b for a,b in zip(xyz_coords[ring_atoms[1]], xyz_coords[ring_atoms[2]]))
                #    new_geo = automol.geom.rotate(new_geo, rotate_xyzs, rotate_angle, xyz_coords[ring_atoms[0]], rotate_hyds)
                #    print(automol.geom.string(new_geo) + '\n')
                if not automol.pot.low_repulsion_struct(geo, new_geo):
                    if automol.geom.is_unique(new_geo, unique_geos, check_dct):
                        unique_zmas.append(nzma)
                        unique_geos.append(new_geo)

        for new_geo in unique_geos:
            print(automol.geom.xyz_string(new_geo))
#        ring_atoms = list(ring_atoms)
#        print(ring_value_dct)
#        defined_atoms = []
#        for da in ring_value_dct:
#            defined_atoms.append(list(coos[da])[0][0])
#
#        plane_idxs = [atom for atom in ring_atoms if atom not in defined_atoms]
#        print(plane_idxs, ring_atoms)
#        alpha_to_plane_idxs = []
#        for atm in plane_idxs:
#            idx = ring_atoms.index(atm)
#            pos = ring_atoms[0 if idx+1 == len(ring_atoms) else idx+1]
#            neg = ring_atoms[idx-1]
#            if pos not in plane_idxs:
#                alpha_to_plane_idxs.append(pos)
#            if neg not in plane_idxs:
#                alpha_to_plane_idxs.append(neg)
#        print('adj to plane', alpha_to_plane_idxs)
#        alpha_to_plane_da = []
#        for da in ring_value_dct:
#            if list(coos[da])[0][0] in alpha_to_plane_idxs:
#                alpha_to_plane_da.append(da)
#        print('adj to plane angles', alpha_to_plane_da)
#
#        beta_to_plane_da = []
#        for atm in alpha_to_plane_idxs:
#            idx = ring_atoms.index(atm)
#            pos = ring_atoms[0 if idx+1 == len(ring_atoms) else idx+1]
#            neg = ring_atoms[idx-1]
#            if pos not in plane_idxs + alpha_to_plane_idxs:
#                for da in ring_value_dct:
#                    if list(coos[da])[0][0] == pos and da not in beta_to_plane_da:
#                        beta_to_plane_da.append(da)
#            if neg not in plane_idxs + alpha_to_plane_idxs:
#                for da in ring_value_dct:
#                    if list(coos[da])[0][0] == neg and da not in beta_to_plane_da:
#                        beta_to_plane_da.append(da)
#        
#        print('beta to plane angles', beta_to_plane_da)
#        val_dct_perms = [{}]
#        if len(beta_to_plane_da) == 0:
#            if len(alpha_to_plane_da) == 1:
#                alpha1 = alpha_to_plane_da[0]
#                q = ring_value_dct[alpha1]
#                val_dct_perms = [{alpha1: q}, {alpha1: -q}]
#            if len(alpha_to_plane_da) == 2:
#                alpha1 = alpha_to_plane_da[0]
#                alpha2 = alpha_to_plane_da[1]
#                q = abs(ring_value_dct[alpha1])
#                print(q)
#                val_dct_perms = [
#                    {alpha1: q, alpha2: -q}, # envelope
#                    {alpha1: -q, alpha2: q}, # envelope (inverted)
#                    {alpha1: q + numpy.pi/8, alpha2: -q}, # envelope2 
#                    {alpha1: -q - numpy.pi/8, alpha2: q}, # envelope2 (inverted)
#                    {alpha2: q + numpy.pi/8, alpha1: -q}, # envelope3
#                    {alpha2: -q - numpy.pi/8, alpha1: q}, # envelope3 (inverted)
#                    {alpha2: q + numpy.pi/16, alpha1: -q + numpy.pi/16}, # half-chair?
#                    {alpha2: - q + numpy.pi/16, alpha1: q + numpy.pi/16}, # half-chair? (inverted)
#                    {alpha2:  - q, alpha1: - q},
#                    {alpha1: -q - numpy.pi/4, alpha2: q}] # half-chair?
#            # if len(beta_to_plane_da) == 1:
#                    
#        elif len(beta_to_plane_da) == 1:
#            for val_dct_perm in val_dct_perms:         
#                if len(beta_to_plane_da) == 1:
#                    alpha1 = alpha_to_plane_da[0]
#                    alpha2 = alpha_to_plane_da[1]
#                    beta1 = beta_to_plane_da[0]
#                    q = abs(ring_value_dct[alpha1])
#                    val_dct_perms = [
#                        {alpha1: q, alpha2: q, beta1: -q},  # III
#                        {alpha1: -q, alpha2: -q, beta1: q}, # III (inverted)
#                        {alpha1: q, alpha2: -q, beta1: 0},  # IV 
#                        {alpha1: -q, alpha2: q, beta1: 0},  # IV (inverted)
#                        {alpha1: q, alpha2: 0, beta1: 0},   # VII
#                        {alpha1: -q, alpha2: 0, beta1: 0},  # VII (inverted)
#                        {alpha1: 0, alpha2: 0, beta1: q},   # VII (inverted2)
#                        {alpha1: 0, alpha2: 0, beta1: -q},  # VII (inverted3)
#                        {alpha1: -q, alpha2: 0, beta1: q},  # IV (inverted2)
#                        {alpha1: -q, alpha2: 0, beta1: q},  # IV (inverted3)
#                        {alpha1: 0, alpha2: -q*2, beta1: q},  # 
#                        {alpha1: q, alpha2: q, beta1: -q*2}  # 
#                    ]
                    #    {alpha1: -q, alpha2: q},
                    #    {alpha1: -q, alpha2: q},
                    #    {alpha1: q + numpy.pi/8, alpha2: -q},
                    #    {alpha1: -q - numpy.pi/8, alpha2: q},
                    #    {alpha2: q + numpy.pi/8, alpha1: -q},
                    #    {alpha2: -q - numpy.pi/8, alpha1: q}]
        # for alpha_da in alpha_to_plane_da:
        #     val_dct_perms_loop = val_dct_perms.copy()
        #     for val_dct_perm in val_dct_perms_loop:
        #         val_dct_perm[alpha_da] = q    
        #         val_dct_perm_j = val_dct_perm.copy()
        #         val_dct_perm_j[alpha_da] = - q
        #         val_dct_perms.append(val_dct_perm_j)
        # for beta_da in beta_to_plane_da:
        #     val_dct_perms_loop = val_dct_perms.copy()
        #     for val_dct_perm in val_dct_perms_loop:
        #         print('heere', val_dct_perm)
        #         val_dct_perm[beta_da] = q
        #         val_dct_perm_j = val_dct_perm.copy()
        #         val_dct_perm_j[beta_da] = - q
        #         val_dct_perm_k = val_dct_perm.copy()
        #         val_dct_perm_k[beta_da] = 0
        #         val_dct_perms.extend([val_dct_perm_j, val_dct_perm_k])
###        print(val_dct_perms)
    
###        for val_dct in val_dct_perms:
###            print(val_dct)
###            new_zma = automol.zmat.set_values_by_name(zma, val_dct, degree=False)
###            new_geo = automol.zmat.geometry(new_zma)
###            print(automol.geom.string(new_geo))
 
        #for i in range(len(ring_atoms)):
    #    subset = ring_atoms[i: min(i + 4, len(ring_atoms))] + ring_atoms[0: max(0, i + 4 - len(ring_atoms))]
    #    print(subset)
    #    found = False
    #    for da in ring_value_dct:
    #        print(coos[da])
    #        if len(list(set(list(coos[da])[0]) & set(subset))) == 3:
    #            found = True
    #    if not found:
    #        print('found plane', subset)
    #        plane_idxs = subset
###    print(plane_idxs)
###    print(ring_value_dct)
    # print('values', ring_value_dct.values())
    # for i, dai in enumerate(ring_value_dct):
    #     for j, daj in enumerate(ring_value_dct):
    #         if i > j:
    #             da_vali = ring_value_dct[dai]
    #             da_valj = ring_value_dct[daj]
    #             print(dai, daj)
    #             new_zma = automol.zmat.set_values_by_name(zma, {dai: 0, daj: 0}, degree=False)
    #             new_geo = automol.zmat.geometry(new_zma)
    #             print(automol.geom.string(new_geo))


if __name__ == '__main__':
    # test__rings()
    # test__ring_systems()
    # test__ring_systems_decomposed_atom_keys()
    test__ring_puckering()
