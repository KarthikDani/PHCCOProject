from setuptools import setup, find_packages

VERSION = '0.1.9'
DESCRIPTION = "Genome Scale Metabolic Modeling - [Currently] Pipeline for building and analyzing context-specific metabolic models"

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='gsmm',
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Karthik Dani',
    author_email='karthikdani14@gmail.com',
    url='https://github.com/KarthikDani/PHCCOProject/tree/main/gsmm',
    packages=find_packages(),
    license='MIT',
    install_requires=[
        'numpy',
        'pandas',
        'seaborn',
        'matplotlib',
        'cobra',
        'corda',
        'scipy',
        'scikit-learn',
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    python_requires='>=3.7'
    )
