""" Install autochem
"""
from distutils.core import setup


setup(name='autochem',
      version='0.1.0',
      packages=['automol',
                'automol.cart',
                'automol.convert',
                'automol.create',
                'automol.dict_',
                'automol.graph',
                'automol.zmatrix',
                'automol.formula',
                'automol.mult',
                'automol.tests',
                'autoread',
                'autoread.zmatrix',
                'autowrite',
                'transformation'],
      package_dir={'automol': 'automol',
                   'autoread': 'autoread',
                   'autowrite': 'autowrite',
                   'transformations', 'transformations'},
      package_data={'automol': ['tests/data/*.txt']})
