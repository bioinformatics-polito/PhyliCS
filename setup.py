from setuptools import find_packages, setup

setup(
    name='phylics',
    packages=find_packages(include=['phylics']),
    version='1.0.0',
    description='Single-cell CNV data analysis toolkit',
    author='Marilisa Montemurro',
    license='AGPL3',
)