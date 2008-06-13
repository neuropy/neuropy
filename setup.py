"""neuropy installation script

to create source distribution and force tar.gz file:
>>> python setup.py sdist --formats=gztar
to create binary distribution:
>>> python setup.py bdist_wininst
"""

from distutils.core import setup
import os

setup(name='neuropy',
      version='0.1',
      license='BSD',
      description='Neuronal spike data and stimulus analysis in Python',
      author='Martin Spacek',
      author_email='mspacek at interchange ubc ca',
      url='http://swindale.ecc.ubc.ca/neuropy',
      #long_description='',
      #packages=['neuropy', 'neuropy.scripts'])
      packages=['neuropy'])
