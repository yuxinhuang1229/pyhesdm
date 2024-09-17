from setuptools import setup
from setuptools.command.install import install
import urllib.request
import os

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        url = 'https://ned.ipac.caltech.edu/NED::LVS/fits/Current/'
        destination = os.path.dirname(__file__)
        if not os.path.exists(destination):
            print(f'Downloading NEDLVS Catalog from {url}...')
            urllib.request.urlretrieve(url, destination)
            print('Download complete.')

setup(
    name='pyhesdm',
    version='0.1',
    packages=['pyhesdm'],
    cmdclass={
        'install': PostInstallCommand,
    }
)