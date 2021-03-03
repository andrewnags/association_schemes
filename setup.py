from setuptools import setup, find_packages


setup(
    name='association_schemes',
    version='0.1',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[],
    tests_require=['pytest', 'pytest-cov', 'flake8', 'mypy'],
    extras_require={
        'tests': ['pytest', 'pytest-cov', 'flake8', 'mypy'],
    },
)
