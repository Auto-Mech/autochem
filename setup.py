""" Install autochem
"""
from distutils.core import setup


setup(name="autochem",
      version="0.2.1",
      packages=["automol",
                "automol.cart",
                "automol.convert",
                "automol.create",
                "automol.dict_",
                "automol.graph",
                "automol.tests",
                "autorxn",
                "autorxn.tests",
                "autoread",
                "autoread.zmatrix",
                "autowrite"],
      package_dir={'automol': "automol", 'autorxn': "autorxn"},
      package_data={'automol': ["tests/data/*.txt"],
                    'autorxn': ["tests/data/*.txt"]},)
