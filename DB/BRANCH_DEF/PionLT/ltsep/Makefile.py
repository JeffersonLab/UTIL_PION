#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2022-06-14 17:54:11 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
'''
from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Set Cuts',
    ext_modules = cythonize("setcut.pyx"),
    zip_safe=False,
)
'''
from distutils.core import Extension, setup
from Cython.Build import cythonize

# define an extension that will be cythonized and compiled
ext = Extension(name="setcut", sources=["setcut.pyx"])
setup(ext_modules=cythonize(ext))
