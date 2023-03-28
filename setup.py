from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

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
    long_description=LONG_DESCRIPTION,
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