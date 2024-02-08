
#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path


this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'src', 'PySulfSat', '_version.py'), encoding='utf-8') as f:
    exec(f.read())

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name="PySulfSat",
    version=__version__,
    author="Penny Wieser",
    author_email="penny.wieser@gmail.com",
    description="PySulfSat",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PennyWieser/PySulfSat",
    package_dir={'': 'src'},  # Optional
    packages=find_packages(where='src'),  # Required

    package_data={
        # Include all pickle files
        "": ["*.pkl"],
    },
    install_requires=[
            'pandas',
            'numpy<2',
            'python-ternary',
            'matplotlib',
            'scikit-learn',
            'scipy',
            'statsmodels',
            'openpyxl',
            ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
