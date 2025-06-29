from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import numpy as np
import os

# Get the path to the pluginsteplight directory
pluginsteplight_dir = os.path.join(os.path.dirname(__file__), '..', 'pluginsteplight')

# Define the Cython extension
extensions = [
    Extension(
        "PythonStep.stl_grid3d_wrapper",
        sources=["PythonStep/stl_grid3d_wrapper.pyx"],
        include_dirs=[
            np.get_include(),
        ],
        extra_compile_args=[
            "-std=c++11",
            "-fPIC",
            "-O3",
        ],
        extra_link_args=[
            "-std=c++11",
        ],
        language="c++",
    )
]

setup(
    name="PythonStep",
    version="0.1.0",
    description="Python wrapper for STL_Grid3D using Cython",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages(),
    ext_modules=cythonize(extensions, compiler_directives={
        'language_level': 3,
        'boundscheck': False,
        'wraparound': False,
        'initializedcheck': False,
        'nonecheck': False,
    }),
    install_requires=[
        "numpy>=1.19.0",
        "cython>=0.29.0",
    ],
    python_requires=">=3.7",
    zip_safe=False,
) 