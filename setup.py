from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="colormycells",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for creating perceptually uniform colormaps for cell types",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/colormycells",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scikit-learn",
        "scipy",
        "colour-science",
    ],
)