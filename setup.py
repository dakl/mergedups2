try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
    
# For making things look nice on pypi:
# try:
#     import pypandoc
#     long_description = pypandoc.convert('README.md', 'rst')
# except (IOError, ImportError):

long_description = 'read merged for PCR duplicates'

setup(name='mergedups2',
    version='0.0.2',
    description='Reduce sequencing error rates by merging PCR duplicates',
    author = 'Daniel Klevebring',
    author_email = 'daniel.klevebring@gmail.com',
    url = 'http://github.com/dakl/mergedups2',
    license = 'MIT License',
    install_requires=['pysam'],
    packages = [
        'mergedups'
    ],
    keywords = [
        'bioinformatics', 
        'pcr',
        'error-correction'
    ],
    scripts = [
        'scripts/mergedups2'
    ],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 3 - Alpha",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    long_description = long_description,
)
