from setuptools import setup, find_packages, Extension
import numpy

ext = Extension(name='dftbephy.extensions', sources=['dftbephy/extensions.pyx'],
                language='c++',
                include_dirs=[numpy.get_include()])

setup(
    name='dftbephy',
    version='0.0.1',
    packages=find_packages(),
    ext_modules=[ext],
    install_requires=['numpy', 'scipy', 'phonopy']
)
