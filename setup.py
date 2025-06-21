from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='bowshockpy',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "astropy",
    ],
    entry_points={
        "console_scripts": [
            "bowshockpy=bowshockpy.genbow:main",
        ],
    },
    long_description=long_description,
    long_descripton_content_type="text/markdown",
)
