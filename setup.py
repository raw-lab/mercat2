import os
import re
import setuptools

main_script = 'bin/mercat2.py'
version = re.search(r'__version__\s+= "(.+)"', open(main_script).read()).group(1)
author = re.search(r'__author__\s+= "(.+)"', open(main_script).read()).group(1)


# recursively load package files
def package_files(directory):
    paths = []
    for (path, _, filenames) in os.walk(directory):
        for filename in filenames:
            if not filename.endswith('.py'):
                paths.append(os.path.join('..', path, filename))
    return paths


setuptools.setup(
    name = "mercat2",
    version = version,
    author = author,
    author_email = "jlfiguer@charlotte.edu",
    description = "versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis",
    long_description = open("README.md", "r").read(),
    long_description_content_type = "text/markdown",
    url = "https://github.com/raw-lab/mercat2",
    scripts = [main_script], # scripts to install to 'bin' path
    packages = ['mercat2_lib'], # list of packages, installed to 'site-packages' folder
    package_dir = dict(mercat2_lib='lib'), # dict with 'package'='relative dir'
    package_data = dict(mercat2_lib=package_files('lib/')), # add non-python data to package, relative paths
    license = "BSD License", # metadata
    platforms = ['Unix'], # metadata
    classifiers = [ # This is the new updated way for metadata, but old way seems to still be used in some of the output
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Unix",
    ],
    python_requires = '>=3.9',
    install_requires = [
        'setuptools',
        'grpcio ==1.43',
        'ray',
        'configargparse',
        'pandas',
        'numpy',
        'humanize',
        'plotly',
        'psutil',
        'dominate',
        'scikit-learn',
        'scikit-bio',
        'scipy',
        'metaomestats',
        'kaleido',

# ray-core ray-dashboard ray-default ray-tune
    ]
)
