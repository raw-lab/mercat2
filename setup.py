import os
import setuptools


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
    version = "0.2",
    author = "Jose Luis Figueroa III, Richard White III",
    author_email = "jlfiguer@uncc.edu",
    description = "versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis",
    long_description = open("README.md", "r").read(),
    long_description_content_type = "text/markdown",
    url = "https://github.com/raw-lab/mercat2",
    scripts = ['bin/mercat2.py'], # scripts to install to 'bin' path
    packages = ['mercat2_lib'], # list of packages, installed to 'site-packages' folder
    package_dir = dict(mercat2_lib='lib'), # dict with 'package'='relative dir'
    package_data = dict(mercat2_lib=package_files('lib/')), # add non-python data to package, relative paths
    license = "BSD License", # metadata
    platforms = ['Unix'], # metadata
    classifiers = [ # This is the new updated way for metadata, but old way seems to still be used in some of the output
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Unix",
    ],
    python_requires = '<3.10',
    install_requires = [
        'setuptools',
        'ray',
        'configargparse',
        'pandas',
        'numpy',
        'humanize',
        'plotly',
        'psutil',
        'dominate',
        'scikit-learn',
    ]
)
