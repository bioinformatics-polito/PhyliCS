from setuptools import find_packages, setup

setup(
    name='phylics',
    packages=find_packages(),
    python_requires='>=3.8, <3.9',
    install_requires=['numpy', 
            'scipy', 
            'pandas', 
            'matplotlib', 
            'seaborn',
            'scikit-learn>=0.24',
            'statsmodels',
            'anndata',
            'typing',
            'umap-learn',
            'IPython',
            'hdbscan',
            'joblib'],
    version='1.0.0',
    description='Single-cell CNV data analysis toolkit',
    author='Marilisa Montemurro',
    license='AGPL3',
)