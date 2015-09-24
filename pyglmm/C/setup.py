# setup.py
from distutils.core import setup, Extension

setup(name="example",
      version="1.0",
      py_modules =['example.py'],
      ext_modules = [
          Extension("_example",
                    ["pyexample.c","example.c"])
      ]
)