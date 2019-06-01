import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AtomicContributions-JaGeo",
    version="1.3",
    author="Janine George",
    author_email="janine.george@uclouvain.be",
    description="Package to display atomic contributions to modes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JaGeo/AtomicContributions",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: BSD 3-Clause License ",
        "Operating System :: OS Independent",
    ],
)
