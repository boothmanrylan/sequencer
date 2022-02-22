from setuptools import setup

setup(
    name='sequencer',
    version='0.0.0',
    url='https://github.com/boothmanrylan/sequencer/',
    description='Search for sequences of events in GEE images.',
    author='Rylan Boothman',
    author_email='boothmanrylan@gmail.com',
    packages=['sequencer'],
    install_requires=['earthengine-api', 'IPython', 'folium']
)
