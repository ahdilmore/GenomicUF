import ast
import re
from setuptools import setup

_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("guf/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

with open("README.md") as f:
    long_description = f.read()

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Software Engineering
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: Only
    Operating System :: Linux
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

classifiers = [s.strip() for s in classes.split("\n") if s]

description = (
    "Constructs per-gene phylogenetic trees from [meta]genomic "
    "sequencing data enabling calculation of UniFrac distance "
    "across multiple gene histories."
)

#standalone = [""]

setup(
    name="GenomicUniFrac", 
    author="Amanda Hazel Dilmore",
    author_email="adilmore@ucsd.edu",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ahdilmore/GenomicUF",
    version=version,
    license="BSD-3-Clause",
    packages=["guf/src"],
    install_requires=[
        "pandas>=1.0.0",
        "numpy",
        "pybedtools",
        "biopython",
        "scikit-bio > 0.5.3",
	"biom-format",
        "unifrac"
    ],
    include_package_data=True,
    package_data={
        "GenomicUniFrac": ["tests"]
    },
    #entry_points={"console_scripts": standalone},
    classifiers=classifiers,
)
