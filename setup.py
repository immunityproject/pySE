"""setup.py -- setup script

Setuptools config
"""

from setuptools import setup

setup(
    name='pyse',
    version='0.0.1',
    packages=['pyse'],
    install_requires=[
        'bigfloat',
        'click'
    ],
    entry_points={
       'console_scripts': [
           'compute_shannon_entropy = pyse.compute_shannon_entropy:main'
       ]
    }
)
