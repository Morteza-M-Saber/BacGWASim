from setuptools import setup

from BacGWASim import __version__, _program

setup(
    name=_program,
    version=__version__,
    author="Masih M. Saber",
    author_email="morteza.mahmoudisaber@gmail.com",
    description="BacGWASim command line interface",
    license="MIT",
    packages=["BacGWASim"],
    entry_points="""
    [console_scripts]
    {program} = BacGWASim.main:main
    """.format(program=_program),
    include_package_data=True,
)
