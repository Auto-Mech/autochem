""" Install autochem
"""
from distutils.core import setup


setup(name="autochem",
      version="0.1.7",
      packages=["automol",
                "automol.constructors",
                "automol.writers",
                "automol.cart",
                "automol.inchi",
                "automol.geom",
                "automol.zmatrix",
                "automol.graph",
                "automol.graph._dict",
                "automol.graph._inchi",
                "automol.graph._stereo",
                "automol.tests",
                "autorxn"],
      package_dir={'automol': "automol"},
      package_data={'automol': ["tests/data/*.txt"]},)
