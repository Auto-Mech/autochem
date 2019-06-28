""" species samples script
"""
import os
import pandas
import automol


if __name__ == '__main__':
    PATH = os.path.dirname(os.path.realpath(__file__))
    TAB = pandas.read_csv(os.path.join(PATH, 'data', 'species.csv'))

    for IDX, ICH in enumerate(TAB['inchi']):
        print(IDX, ICH)
        CGR = automol.inchi.graph(ICH)
        for SGR in automol.graph.stereomers(CGR):
            SICH = automol.graph.inchi(SGR)
            print(SICH)
