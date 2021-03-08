""" Install automol
"""
from distutils.core import setup


setup(name='automol',
      version='0.1.1',
      packages=['automol',
                'automol.cart',
                'automol.convert',
                'automol.create',
                'automol.dict_',
                'automol.embed',
                'automol.etrans',
                'automol.formula',
                'automol.geom',
                'automol.graph',
                'automol.instab',
                'automol.intmol',
                'automol.mult',
                'automol.pot',
                'automol.prop',
                'automol.reac',
                'automol.rotor',
                'automol.zmat',
                'automol.tests',
                'automol.util',
                'automol.util.dict_',
                'phydat',
                'transformations'],
      package_dir={'automol': 'automol',
                   'phydat': 'phydat',
                   'transformations': 'transformations'},
      package_data={'automol': ['tests/data/*.txt']})
