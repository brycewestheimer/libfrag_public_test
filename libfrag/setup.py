from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import os

# Read the version from the C++ header
version_file = os.path.join("include", "libfrag", "libfrag_version_major.hpp")
version_major = 0
version_minor = 1
version_patch = 0

# Try to read version from headers if available
try:
    with open(os.path.join("include", "libfrag", "libfrag_version_major.hpp"), 'r') as f:
        for line in f:
            if "#define LIBFRAG_VERSION_MAJOR" in line:
                version_major = int(line.split()[-1])
    with open(os.path.join("include", "libfrag", "libfrag_version_minor.hpp"), 'r') as f:
        for line in f:
            if "#define LIBFRAG_VERSION_MINOR" in line:
                version_minor = int(line.split()[-1])
    with open(os.path.join("include", "libfrag", "libfrag_version_patch.hpp"), 'r') as f:
        for line in f:
            if "#define LIBFRAG_VERSION_PATCH" in line:
                version_patch = int(line.split()[-1])
except FileNotFoundError:
    pass

__version__ = f"{version_major}.{version_minor}.{version_patch}"

# Read long description from README
with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="libfrag",
    version=__version__,
    author="Your Name",  # Update this
    author_email="your.email@example.com",  # Update this
    description="A quantum chemistry library for fragment-based calculations",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/yourusername/libfrag",  # Update this
    packages=find_packages(where="python/module"),
    package_dir={"": "python/module"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",  # Update if different
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.19.0",
    ],
    extras_require={
        "dev": ["pytest", "pybind11>=2.6.0"],
    },
)
