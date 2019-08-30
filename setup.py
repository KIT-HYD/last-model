from setuptools import setup, find_packages


def requirements():
    with open('requirements.txt') as f:
        return f.read().split('\n')


def version():
    with open('VERSION') as f:
        return f.read().strip()


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name="last_model",
    author="Mirko Mälicke, Alexander Sternagel",
    author_email="mirko.maelicke@kit.edu",
    install_requires=requirements(),
    version=version(),
    description="LAST model Python implementation",
    long_description=readme(),
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False
)