# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:59:20 2016

@author: Antoine
"""

from distutils.core import setup

setup(
    name='Solar_mod',
    version='0.1.0',
    author='Antoine',
    author_email='antoine.letarte@gmail.com',
    packages=['solar_mod'],
    scripts=['bin/test.py'],
#    url='http://pypi.python.org/pypi/Solar/',
#    license='LICENSE.txt',
    description='ENR835',
    long_description=open('README.txt').read(),
    install_requires=[
        "Django >= 1.1.1",
        "caldav == 0.1.4",
    ],
)


