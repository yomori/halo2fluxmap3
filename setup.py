from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension("halo2fluxmap3.flux2map"
              ,["halo2fluxmap3/flux2map.pyx"],
              include_dirs=[np.get_include()]
    ),]

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(parent_package,top_path)
    return config

if __name__ == "__main__":
      from numpy.distutils.core import setup
      setup(name='halo2fluxmap3',
      configuration=configuration,
      ext_modules=cythonize(ext_modules),
      include_dirs=[np.get_include()],
      version='0.1',
      description='for converting halo catalogs to maps for a given hod',
      url='http://github.com/marcelo-alvarez/halo2fluxmap',
      author='Marcelo Alvarez',
      author_email='marcelo.alvarez@berkeley.edu',
      license='GNU General Public License v3.0',
      packages=['halo2fluxmap3'],
      zip_safe=False)

