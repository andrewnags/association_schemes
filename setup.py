from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy as np
from sage.env import sage_include_directories


setup(
    name='association_schemes',
    version='0.1',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=["numpy"],
    setup_requires=[
        "setuptools",
        "Cython",
    ],
    tests_require=['pytest', 'pytest-cov', 'flake8', 'mypy'],
    extras_require={
        'tests': ['pytest', 'pytest-cov', 'flake8', 'mypy'],
    },
    include_dirs=[
        np.get_include(),
    ]  + sage_include_directories(),
    ext_modules=cythonize("src/**/*.pyx", language_level=3),
)
