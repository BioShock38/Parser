from setuptools import setup, find_packages

setup(name='polymnie',
      version='0.1',
      description='Parser and writer in different bioinformatic format',
      author='Thomas Dias-Alves',
      packages=find_packages(),
      test_suite='pyhaplophase.tests',
      zip_safe=False)
