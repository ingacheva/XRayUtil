from setuptools import setup, find_packages
from os.path import join, dirname

setup(name='XRayUtil',
      version='1.0',
      long_description=open(join(dirname(__file__), 'README.md')).read(),
      url='https://github.com/ingacheva/XRayUtil',
      author='Anastasia Ingacheva',
      author_email='ingacheva@gmail.com',
      license='MIT',
      packages=find_packages(),
      dependency_links=['https://anaconda.org/conda-forge/xraylib/files#egg=xraylib-3.3.0'],
      include_package_data=True,
      zip_safe=False)