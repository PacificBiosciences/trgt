import codecs
import os.path
from setuptools import setup
import subprocess

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

VERSION = "0.0.1"

setup(
    name="TRGT",
    version=VERSION,
    author="Egor Dolzenhko",
    author_email="",
    url="https://github.com/PacificBiosciences/trgt",
    packages=["trgt", "trgt/database"],
    license="BSD 3-Clause Clear License",
    description="Tandem repeat genotyping and visualization from PacBio HiFi data",
    include_package_data=True,
    long_description=open("README.md", encoding="UTF-8").read(),
    long_description_content_type="text/markdown",
    entry_points={
      "console_scripts": [
         "trgt = trgt.__main__:main"
      ]
    },
    install_requires=[
        "truvari>=4.0",
        "pysam>=0.15.2",
        "pyarrow==11.0.0",
    ],
)
