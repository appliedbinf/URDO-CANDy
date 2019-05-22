#!/usr/bin/env python
try:
  import os
  from setuptools import setup, find_packages
except ImportError:
  from distutils.core import setup
from os import path
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md')) as fh:
      long_description_text = fh.read()


setup(
  name = "CANDY",
  scripts = ['candy.py',],
  version = "0.1.0",
  description = 'Fast k-mer based tool for alignment and assembly-free amplicon classification.',
  long_description=long_description_text,
  long_description_content_type="text/markdown",
  author = 'Anna Gaines',
  author_email = 'anna.gaines@ihrc.com',
  url = 'https://github.com/ar0ch/URDO-CANDy',
  classifiers = [
      'Programming Language :: Python :: 3.6',
  ],
  install_requires=['lxml'],
)
