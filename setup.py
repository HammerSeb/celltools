from gettext import find
from setuptools import setup, find_packages

setup(name='celltools',
      version='0.0.2',
      description='Tools to manipulate atoms and molecules in a unit cell',
      # long_description=open('README.md').read(),
      # long_description_content_type ='text/markdown' ,
      license= 'MIT',
      author='Sebastian Hammer',
      author_email='sebastian.hammer@mail.mcgill.ca',
      url='https://github.com/HammerSeb/celltools',
      packages= find_packages(),
      install_requires=['numpy', 'crystals', 'matplotlib'],
      python_requires='>=3.6',
      include_package_data = True
     )