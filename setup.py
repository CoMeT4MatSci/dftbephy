from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
  name='dftbephy',
  ext_modules=cythonize('dftbephy/extensions.pyx', annotate=True) 
)
