from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "pre-BotC model",
    ext_modules = cythonize('*.pyx'), # accepts a glob pattern
)
