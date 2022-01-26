""" Install autochem
"""

from distutils.core import setup


setup(
    name='autochem',
    version='0.13.1',
    packages=[
        'automol',
        # L1
        'automol.util',
        'automol.util.dict_',
        'automol.mult',
        'automol.formula',
        'automol.prop',
        'automol.embed',
        # L2
        'automol.geom.base',
        'automol.graph.base',
        'automol.zmat.base',
        # L3
        'automol.extern',
        'automol.inchi.base',
        # L4
        'automol.graph',
        'automol.geom',
        'automol.inchi',
        'automol.zmat',
        # L5
        'automol.pot',
        'automol.etrans',
        'automol.reac',
        'automol.rotor',
        # other
        'autoreact',
        'phydat',
        'transformations'
    ],
    package_dir={
        'automol': 'automol',
        'autoreact': 'autoreact',
        'phydat': 'phydat',
        'transformations': 'transformations'
    },
    package_data={
        'automol': ['tests/data/*'],
    }
)
