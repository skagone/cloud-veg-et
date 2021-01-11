from setuptools import setup

setup(name='vegetLib',
      maintainer='Tony Butzer',
      maintainer_email='tonybutzer@gmail.com',
      version='1.0.0',
      description='Classes for et model water balance',
      packages=[
          'vegetLib',
      ],
      install_requires=[
          'boto3',
      ],

)
