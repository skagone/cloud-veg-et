from setuptools import setup

setup(name='playLib',
      maintainer='Tony Butzer',
      maintainer_email='tonybutzer@gmail.com',
      version='1.0.3',
      description='helper functions for et model water balance',
      packages=[
          'playLib',
      ],
      install_requires=[
          'boto3',
      ],

)
