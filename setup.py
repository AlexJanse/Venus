#!/usr/bin/env python
import os
import sys
from setuptools import setup, find_packages

if sys.version_info < (3, 6):
    sys.exit(
        'Python < 3.6 is not supported. You are using Python {}.{}.'.format(
            sys.version_info[0], sys.version_info[1])
    )

with open('requirements.txt', 'r') as f:
    required_packages = f.read().splitlines()

with open(os.path.dirname(os.path.realpath(__file__)) + '/__version__', 'r') as ver:
    version = ver.readline().strip()
setup(
    name="Venus",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'module_detection = src.module_detection.module_detection:main',
        ]
    },
    zip_safe=False,
    python_requires='>=3.6.0',
    install_requires=required_packages,
    version=version
)
