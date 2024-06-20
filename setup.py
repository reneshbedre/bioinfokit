from setuptools import setup
import bioinfokit

with open("README.md", "r", encoding="UTF-8") as fh:
    long_description = fh.read()

setup(name='bioinfokit',
      version=bioinfokit.__version__,
      # metadata
      author='Renesh Bedre',
      author_email='reneshbe@gmail.com',
      maintainer='Renesh Bedre',
      description='Bioinformatics data analysis and visualization toolkit',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT',
      url='https://github.com/reneshbedre/bioinfokit',
      packages=['bioinfokit',],
      zip_safe=False,
      classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
       ],

      install_requires=['pandas', 'numpy', 'matplotlib', 'scipy', 'scikit-learn', 'seaborn', 'matplotlib-venn',
                        'tabulate', 'statsmodels', 'textwrap3', 'adjustText'],
      )
