""" Install autochem
"""
from distutils.core import setup


setup(name="autochem",
      version="0.2.4",
      packages=["automol",
                "automol.cart",
                "automol.convert",
                "automol.create",
                "automol.dict_",
                "automol.graph",
                "automol.tests",
                "autoread",
                "autoread.zmatrix",
                "autowrite"],
      package_dir={'automol': "automol"},
      package_data={'automol': ["tests/data/*.txt"]})
