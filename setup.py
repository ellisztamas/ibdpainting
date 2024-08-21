from setuptools import setup
import codecs
import os.path
from pathlib import Path

# Functions to pull the package version from init.py
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

# read the contents of your README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='ibdpainting',
    version=get_version("ibdpainting/__init__.py"),
    description='Identify parents of a crossed individual by comparing identity in windows across their genomes',
    url='https://github.com/ellisztamas/ibdpainting',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Tom Ellis',
    author_email='thomas.ellis@gmi.oeaw.ac.at',
    license='MIT',
    packages=['ibdpainting'],
    install_requires=[
        'numpy', 'pandas', 'plotly', 'h5py', 'scikit-allel', 'pytest', 'kaleido'
      ],
    zip_safe=False,
    entry_points = {
        'console_scripts': ['ibdpainting=ibdpainting.command_line:main'],
    }
    )

