from setuptools import setup, find_packages

setup(name='polymnie',
      version='0.1',
      description='Parser and writer in different bioinformatic format',
      author='Thomas Dias-Alves',
      packages=find_packages(),
      test_suite='polymnie.tests',
      install_requires=[
          'vcfnp'
      ],
      zip_safe=False)
