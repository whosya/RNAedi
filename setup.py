from setuptools import setup
from rnaediting._version import __version__

# with open("README.md", 'r') as f:
#     long_description = f.read()

# with open('requirements.txt') as f:
#     required = f.read().splitlines()

setup(
    name='RNAedi',
    version=__version__,
    packages=['rnaediting'],    
    url='https://github.com/whosya/rnaediting',
    license='GNU GPL v3.0',
    author='Shuaiya Hu',
    author_email='whosy@stu.njau.edu.cn',
    description='Command line tool for rnaediting analysis.',
    # long_description=long_description,
    # package_data={
    #     'ksrates': ['ks.mplstyle']
    # },
    py_modules=['RNAedi_cli'],
    #install_requires=[required],
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "Operating System :: OS Independent"
    # ],
    python_requires='>=3.8',
    entry_points={
        "console_scripts": [
            "RNAedi = RNAedi_cli:cli"
        ]
    },
)
