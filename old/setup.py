
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'PrimeTest',
  ext_modules=[
    Extension('c_reactionKinetics', ['c_reactionKinetics.pyx'])
    ],
  cmdclass = {'build_ext': build_ext}
)
