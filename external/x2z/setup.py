""" Install pyx2z

(first run `cmake; make install` to generate the pyx2z shared library)
"""
from distutils.core import setup


setup(name="pyx2z",
      version="0.1.3",
      packages=[''],
      package_dir={'': '.'},
      package_data={'': ["pyx2z*.so"]},)
