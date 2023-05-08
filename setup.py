from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.1'
DESCRIPTION = 'Package for the generation of non-polar and stoichiometric surfaces from ionic structures.'
LONG_DESCRIPTION = 'Package for the generation of non-polar and stoichiometric surfaces from ionic structures.'

# Setting up
setup(
    name="polycleaver",
    version=VERSION,
    author="Eric Mates-Torres",
    author_email="<eric.mates@uab.cat>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=['pymatgen'],
    keywords=['python', 'slab', 'surface', 'miller', 'ionic', 'mineral'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
