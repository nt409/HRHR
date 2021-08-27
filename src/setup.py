from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()



setup(name='hrhr',
      version='0.0.3',
      description='HRHR model',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/nt409/HRHR',
      author='Nick Taylor',
      author_email='nt409@cam.ac.uk',
      license='MIT',
      packages=['model'],
    #   packages=find_packages(),
      zip_safe=False)