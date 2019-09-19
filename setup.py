from setuptools import setup, find_packages
from os import path

__version__ = '0.1.9'

setup(
    name='vdj',
    version=__version__,
    description='Python utilities used for analysis of vdj data.',
    license='MIT',
    classifiers=[
        'Development Status :: beta ',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
    ],
    author="Soichi Hirokawa and Griffin Chure",
)