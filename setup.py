from setuptools import setup

from BacGWASim import __version__, _program

setup(
    name=_program,
    version=__version__,
    packages=["BacGWASim"],
    description="BacGWASim command line interface",
    entry_points="""
    [console_scripts]
    {program} = BacGWASim.main:main
    """.format(program=_program),
    include_package_data=True,
)
