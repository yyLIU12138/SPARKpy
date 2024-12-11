from setuptools import Command, find_packages, setup

__lib_name__ = "SPARKpy"
__lib_version__ = "1.0.0"
__description__ = "Python implementation of SPARK following the original paper"
__url__ = "https://github.com/yyLIU12138/SPARKpy"
__author__ = "Yuyao Liu"
__author_email__ = "2645751574@qq.com"

setup(
    name = __lib_name__,
    version = __lib_version__,
    description = __description__,
    url = __url__,
    author = __author__,
    author_email = __author_email__,
    packages = ['SPARKpy'],
    zip_safe = False,
    include_package_data = True,
)
