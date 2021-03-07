from setuptools import find_packages, setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()

setup(
    name='phylics',
    version='1.0.0',
    description='Single-cell CNV data analysis toolkit',
    #long_description = 'README.md',
    url='https://github.com/bioinformatics-polito/PhyliCS',
    author='Marilisa Montemurro',
    maintainer='Elena Grassi, Gianvito Urgese',
    maintainer_email='elena.grassi@irccs, gianvito.urgese@polito.it',
    author_email='marilisa.montemurro@polito.it', 
    classifiers = [
            'License :: OSI Approved :: AGPL3 License',
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Framework :: Jupyter',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='AGPL3',
    packages=find_packages(),
    #package_data={'phylics': ['datasets/*']},
    python_requires='>=3.8, <3.9',
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