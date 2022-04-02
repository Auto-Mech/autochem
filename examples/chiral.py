""" Basic demo of distance geometry functionality
"""
# import sys
import numpy
import automol


def main(ich):
    """ main function
    """
    # 1. Print inchi
    print('New run')
    print(ich)
    print(automol.inchi.smiles(ich))

    geo_in = automol.inchi.geometry(ich)
    gra = automol.geom.graph(geo_in)

    geo = automol.graph.embed.geometry(gra)

    print("Geometry in:")
    print(automol.geom.string(geo_in))
    print()

    print("Geometry:")
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
    automol.smiles.inchi(r'CC[C@@H](C)O'),
    automol.smiles.inchi(r'CC[C@H](C)O'),
    # automol.smiles.inchi(r'C1C[C@@H]2CN[C@H]1OO2'),
    automol.smiles.inchi(r'F[C@@H](Cl)[C@@H](F)Cl'),
    automol.smiles.inchi(r'F[C@@H](Cl)[C@H](F)Cl'),
    automol.smiles.inchi(r'F[C@H](Cl)[C@H](F)Cl'),
    automol.smiles.inchi(r'FC([C@@H](F)Cl)[C@@H](F)Cl'),
    automol.smiles.inchi(r'FC([C@H](F)Cl)[C@H](F)Cl'),
    automol.smiles.inchi(r'F[C@@H](Cl)[C@@H](F)[C@H](F)Cl'),
    automol.smiles.inchi(r'F[C@@H](Cl)[C@H](F)[C@H](F)Cl'),
    automol.smiles.inchi(r'C=C[C@@H](C)[C@H]([O])/C=C/C'),
    automol.smiles.inchi(r'C=C[C@@H](C)[C@H]([O])/C=C\C'),
    automol.smiles.inchi(r'C=C[C@@H](C)[C@@H]([O])/C=C/C'),
    automol.smiles.inchi(r'C=C[C@@H](C)[C@@H]([O])/C=C\C'),
    automol.smiles.inchi(r'C=C[C@H](C)[C@H]([O])/C=C/C'),
    automol.smiles.inchi(r'C=C[C@H](C)[C@H]([O])/C=C\C'),
    automol.smiles.inchi(r'C=C[C@H](C)[C@@H]([O])/C=C/C'),
    automol.smiles.inchi(r'C=C[C@H](C)[C@@H]([O])/C=C\C'),
]

NTOTAL = 0
NFAIL = 0

for ICH in ICHS:
    # for seed in [1]:
    for seed in range(10):
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
