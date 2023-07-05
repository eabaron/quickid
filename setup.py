# -*- coding: utf-8 -*-
import setuptools
from setuptools import setup

setup(name='quickid', version='0.1',
      description='Quick Line IDs for Observation',
      url='https://github.com/eabaron/quickid',
      author='Eddie Baron',
      author_email='ebaron@psi.edu',
      license='GPL-v3',
      packages=setuptools.find_packages(),
      include_package_data=True,
      install_requires=['numpy', 'mendeleev', 'pandas', 'astropy'],
      optional=['matplotlib','pylab'],
      classifiers=['Programming Language :: Python',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Intended Audience :: Science/Research',
                   'Development Status :: 3 - Alpha']
      )
