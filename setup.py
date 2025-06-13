from setuptools import setup, find_packages

setup(
    name="optimization-projects",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "jupyter>=1.0.0",
        "pandas>=1.3.0",
        "cvxpy>=1.1.0",
        "scikit-learn>=0.24.0",
        "pytest>=6.2.0",
    ],
    author="Salar Nouri",
    author_email="your.email@example.com",
    description="A collection of optimization, game theory, and system identification implementations",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/Optimization",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
) 
