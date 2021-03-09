from setuptools import find_packages, setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README_PyPI.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='phylics',
    version='1.0.5',
    description='Single-cell CNV data analysis toolkit',
    #long_description = 'README.md',
    url='https://github.com/bioinformatics-polito/PhyliCS',
    author='Marilisa Montemurro',
    author_email='marilisa.montemurro@polito.it', 
    classifiers = [
            'License :: OSI Approved :: GNU Affero General Public License v3',
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Framework :: Jupyter',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='AGPL3',
    packages=find_packages(),
    #package_data={'phylics': ['datasets/*']},
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires='>=3.7, <3.9',
    install_requires=['numpy>=1.19.5', 
            'scipy>=1.6.0', 
            'pandas>=1.1.3', 
            'matplotlib>=3.3.1', 
            'seaborn>=0.11.1',
            'scikit-learn>=0.24',
            'statsmodels>=0.12.0',
            'anndata>=0.7.5',
            'typing',
            'umap-learn>=0.4.6',
            'IPython>=7.19.0',
            'hdbscan>=0.8.26',
            'joblib>=1.0.0'],
    
)