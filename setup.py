from setuptools import find_packages, setup


def get_version():
    """
    Gets the version from version.py
    """
    version_file = "bowshockpy/version.py"
    with open(version_file, "r", encoding="utf-8") as f:
        exec(compile(f.read(), version_file, "exec"))
    return locals()["__version__"]


with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="bowshockpy",
    version=get_version(),
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "astropy",
        "photutils",
    ],
    entry_points={
        "console_scripts": [
            "bowshockpy=bowshockpy.genbow:main",
        ],
    },
    author="Guillermo Blazquez-Calero",
    maintainer="Guillermo Blazquez-Calero",
    long_description=long_description,
    long_description_content_type="text/markdown",
)
