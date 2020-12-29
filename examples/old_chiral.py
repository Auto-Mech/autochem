""" Basic demo of distance geometry functionality
"""
import sys
import numpy
import automol


def main(ich):
    """ main function
    """
    # 1. Print inchi
    print('New run')
    print(ich)
    print(automol.inchi.smiles(ich))

    # 2. Generate graph and sorted list of atom keys
    geo = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo)
    keys = sorted(automol.graph.atom_keys(gra))
    syms = list(map(automol.graph.atom_symbols(gra).__getitem__, keys))

    # 3. Generate bounds matrices
    lmat, umat = automol.graph.embed.distance_bounds_matrices(gra, keys)
    pla_dct = automol.graph.embed.planarity_constraint_bounds(gra, keys)
    chi_dct = automol.graph.embed.chirality_constraint_bounds(gra, keys)
    print("Lower bounds matrix:")
    print(numpy.round(lmat, 1))
    print("Upper bounds matrix:")
    print(numpy.round(umat, 1))
    print()
    print(chi_dct)

    # 4. Sample a geometry from the bounds matrices
    xmat = automol.embed.sample_raw_distance_coordinates(lmat, umat, dim4=True)
    geo_init = automol.embed.geometry_from_coordinates(xmat, syms)

    # 5. Clean up the sample's coordinates
    conv_ = automol.graph.embed.qualitative_convergence_checker_(gra, keys)
    xmat, _ = automol.embed.cleaned_up_coordinates(
        xmat, lmat, umat, pla_dct=pla_dct, chi_dct=chi_dct, chi_flip=True,
        conv_=conv_)
    geo = automol.embed.geometry_from_coordinates(xmat, syms)

    # 6. Print the largest errors
    dmat = automol.embed.distance_matrix_from_coordinates(xmat)
    err_dct = automol.embed.greatest_distance_errors(dmat, lmat, umat)
    sp_dct = automol.graph.atom_shortest_paths(gra)
    for (key1, key2), err_val in err_dct.items():
        print('\tError:', err_val)
        if key2 not in sp_dct[key1]:
            print('\nNot connected:', key1, key2)
        else:
            print('\tPath:', sp_dct[key1][key2])

    # 7. Print geometries
    print("Sample geometry:")
    print(automol.geom.string(geo_init))
    print()

    print("Cleaned up geometry:")
    print(automol.geom.string(geo))
    print()

    # 8. Check the connectivity
    gra2 = automol.geom.graph(geo)
    print("Is the graph the same?", 'Yes' if gra == gra2 else 'No')
    assert gra == gra2, (
        'Original graph:\n{:s}\nFinal graph:\n{:s}'
        .format(automol.graph.string(gra), automol.graph.string(gra2)))

    print(automol.graph.string(gra))

    ich = automol.geom.inchi(geo)
    smi = automol.inchi.smiles(ich)
    print(ich)
    print(smi)
    print('End run')
    print()


ICHS = [
    # automol.smiles.inchi(r'CC[C@@H](C)O'),
    # automol.smiles.inchi(r'CC[C@H](C)O'),
    # automol.smiles.inchi(r'C1C[C@@H]2CN[C@H]1OO2'),
    # automol.smiles.inchi(r'F[C@@H](Cl)[C@@H](F)Cl'),
    # automol.smiles.inchi(r'F[C@@H](Cl)[C@H](F)Cl'),
    # automol.smiles.inchi(r'F[C@H](Cl)[C@H](F)Cl'),
    # automol.smiles.inchi(r'FC([C@@H](F)Cl)[C@@H](F)Cl'),
    # automol.smiles.inchi(r'FC([C@H](F)Cl)[C@H](F)Cl'),
    # automol.smiles.inchi(r'F[C@@H](Cl)[C@@H](F)[C@H](F)Cl'),
    # automol.smiles.inchi(r'F[C@@H](Cl)[C@H](F)[C@H](F)Cl'),
    automol.smiles.inchi(r'C=C[C@@H](C)[C@H]([O])/C=C/C'),
    # automol.smiles.inchi(r'C=C[C@@H](C)[C@H]([O])/C=C\C'),
    # automol.smiles.inchi(r'C=C[C@@H](C)[C@@H]([O])/C=C/C'),
    # automol.smiles.inchi(r'C=C[C@@H](C)[C@@H]([O])/C=C\C'),
    # automol.smiles.inchi(r'C=C[C@H](C)[C@H]([O])/C=C/C'),
    # automol.smiles.inchi(r'C=C[C@H](C)[C@H]([O])/C=C\C'),
    # automol.smiles.inchi(r'C=C[C@H](C)[C@@H]([O])/C=C/C'),
    # automol.smiles.inchi(r'C=C[C@H](C)[C@@H]([O])/C=C\C'),
]

NTOTAL = 0
NFAIL = 0

for ICH in ICHS:
    # for seed in range(10):
    for seed in [1]:
        NTOTAL += 1
        numpy.random.seed(seed)
        main(ICH)
        # try:
        #     main(ICH)
        # except Exception as err:
        #     print(err)
        #     print("Seed: {:d}".format(seed))
        #     SMI = automol.inchi.smiles(ICH)
        #     print("Smiles: {:s}".format(SMI))
        #     NFAIL += 1
        #     sys.exit()

print()
print("Total:", NTOTAL)
print("Failed:", NFAIL)
print("%:", NFAIL/NTOTAL*100.)
