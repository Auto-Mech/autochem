""" Basic demo of distance geometry functionality
"""
# import sys
import numpy
import automol


def main():
    """ main function
    """
    # 1. Choose molecule
    # ich = automol.smiles.inchi('CC[C@@H](C)O')
    # ich = automol.smiles.inchi('CC[C@H](C)O')
    ich = automol.smiles.inchi('C1C[C@@H]2CN[C@H]1OO2')

    # 2. Generate graph and sorted list of atom keys
    geo = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo)
    keys = sorted(automol.graph.atom_keys(gra))
    syms = list(map(automol.graph.atom_symbols(gra).__getitem__, keys))

    # 3. Generate bounds matrices
    lmat, umat = automol.graph.embed.distance_bounds_matrices(gra, keys)
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
    xmat, _ = automol.embed.cleaned_up_coordinates(
        xmat, lmat, umat, chi_dct=chi_dct, chi_flip=True)
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

    print(automol.graph.string(gra))

    ich = automol.geom.inchi(geo)
    smi = automol.inchi.smiles(ich)
    print(ich)
    print(smi)


for ntry in range(1):
    numpy.random.seed(ntry)
    try:
        main()
    except Exception as err:
        print(err)
        print("Seed: {:d}".format(ntry))
        break
