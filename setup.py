from setuptools import setup, find_packages

setup(name='celltools',
      version='0.3.0',
      description='Tools to manipulate atoms and molecules in a unit cell',
      long_description=open('README.md').read(),
      long_description_content_type ='text/markdown',
      license= 'MIT',
      author='Sebastian Hammer',
      author_email='sebastian.hammer@mail.mcgill.ca',
      url='https://github.com/HammerSeb/celltools',
      packages= find_packages(),
      install_requires=['numpy', 'scipy', 'pyqtgraph', 'crystals', 'scikit-ued'],
      python_requires='>=3.9',
     )