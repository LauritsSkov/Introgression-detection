from setuptools import setup

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='hmmix', 
    version = '0.7.3',
    description='Find introgressed segments',
    py_modules=['bcf_vcf', 'helper_functions', 'hmm_functions', 'main', 'make_mutationrate', 'make_test_data'],
    package_dir={'': 'src'},
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: MIT License',
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = 'https://github.com/LauritsSkov/Introgression-detection',
    author = 'Laurits Skov and Moises Coll Macia',
    author_email='lauritsskov2@gmail.com',
    entry_points = {
    'console_scripts': [
        'hmmix = main:main'
    ]},
    install_requires=[
          'numpy>=1.15',
          'numba'
      ],   
)

