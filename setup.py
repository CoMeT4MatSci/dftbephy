from setuptools import setup, find_packages, Extension
import numpy

setup_args = dict(
    packages=find_packages(
        where='dftbephy',  # '.' by default
    ),
    ext_modules = [Extension(name='dftbephy.extensions', sources=['dftbephy/extensions.pyx'],
                    language='c++',
                    include_dirs=[numpy.get_include()])]
    )

setup(**setup_args)

