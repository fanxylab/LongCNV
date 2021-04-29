#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# cython: language_level=3
'''
***********************************************************
* @File    : setup.py
* @Author  : Zhou Wei                                     *
* @Date    : 2020/12/11 20:42:39                          *
* @E-mail  : welljoea@gmail.com                           *
* @Version : --                                           *
* You are using the program scripted by Zhou Wei.         *
* If you find some bugs, please                           *
* Please let me know and acknowledge in your publication. *
* Thank you!                                              *
* Best wishes!                                            *
***********************************************************
'''

import sys
from os.path import dirname, join
from glob import glob
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy


extensions = [
    Extension("EcMagiccube", 
        sources=['Script/EcMagiccube.pyx'], 
        language='c',
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        include_dirs=[numpy.get_include()],
        #include_path = [numpy.get_include()],
        library_dirs=[],
        libraries=[],
        extra_compile_args=[],
        extra_link_args=[]),
]

'''
setup(
    name="LongCNV",
    ext_modules=cythonize(extensions, compiler_directives={'language_level' : "3"}),
)
#python setup.py build_ext --inplace
'''

setup_args = {}
# Dependencies for easy_install and pip:
install_requires=[
        'joblib >= 0.13.2',
        'matplotlib >= 3.0.3',
        'numpy >= 1.16.4',
        'pandas >= 0.24.2',
        'Cython >= 0.29.21',
        'numba >= 0.50.1',
]

DIR = (dirname(__file__) or '.')

setup_args.update(
    name='LongCNV',
    version='0.1',
    description=__doc__,
    author='Wei Zhou',
    author_email='welljoea@gmail.com',
    maintainer='Wei Zhou',
    #license = "MIT Licence", 
    url='https://github.com/WellJoea/LongCNV',

    packages = find_packages(where='.', exclude=(), include=('*',)),
    include_package_data = True,

    platforms = "any",
    scripts=[join(DIR, 'LongCNV.py')] + glob(join(DIR, 'Scripts/*.py')),
    ext_modules=cythonize(extensions, compiler_directives={'language_level' : "3"}), 
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
    ],

)

setup(**setup_args)
