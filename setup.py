from setuptools import setup
from setuptools.command.install import install
import urllib.request
import os

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        url = 'https://ned.ipac.caltech.edu/NED::LVS/fits/Current/'
        destination = f"{os.path.dirname(__file__)}/pyhesdm/NEDLVS_20210922_v2.fits"
        if not os.path.exists(destination):
            print(f'Downloading NEDLVS Catalog from {url}...')
            urllib.request.urlretrieve(url, destination)
            print('Download complete.')

setup(
    name='pyhesdm',
    version='0.1',
    description='Local Universe Dispersion Measure Model Computed from HESTIA Simulation',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Yuxin Huang, Khee-Gan Lee, Noam Libeskind, Sunil Simha, Aurélien Valade and J. X. Prochaska',
    author_email='mochafhxy@gmail.com',
    url='https://github.com/yuxinhuang1229/pyhesdm',
    packages=['pyhesdm'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    cmdclass={
        'install': PostInstallCommand,
    }
)