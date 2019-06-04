import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="AtomicContributions_JaGeo",
    version="1.5.0",
    author="Janine George",
    author_email="janine.george@uclouvain.be",
    description="Package to display atomic contributions to modes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JaGeo/AtomicContributions",
    packages=setuptools.find_packages(),
    setup_requires=['setuptools>=18.0'],
    install_requires=["numpy>=1.14.3","phonopy>=2.1.1","matplotlib>=3.0.2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)
