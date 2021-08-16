from setuptools import setup, Extension, find_packages
import os
import numpy
import glob
import sys
from shutil import copyfile

srcFiles = [
    "/Users/csayres/installTCC/coordConv/python/coordConv/coordConvLib.i",
]

# slalibFiles = glob.glob("/Users/csayres/installTCC/slalib_2013-04-25/src/*.c")

# srcFiles += slalibFiles

convFiles = glob.glob("/Users/csayres/installTCC/coordConv/src/*.cc")
srcFiles += convFiles

includeList = [
    "/Users/csayres/installTCC/coordConv/python",
    "Users/csayres/installTCC/coordConv/include",
    numpy.get_include(),
    "/Users/csayres/installTCC/ndarray/8.0.0.0+1/python",
    "/Users/csayres/installTCC/ndarray/8.0.0.0+1/include",
    "/Users/csayres/installTCC/coordConv/include",
    "/Users/csayres/condatcc/include/python2.7",
    "/Users/csayres/installTCC/eigen/3.1.1+2/include",
    "/Users/csayres/installTCC/boost_1_55_0",
    "/Users/csayres/installTCC/slalib_2013-04-25/install/include",
]

extra_compile_args = ["--std=c++11", "-fPIC", "-v", "-O3"]
extra_link_args = None
if sys.platform == 'darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.9']
    extra_link_args = ["-v", '-mmacosx-version-min=10.9']

swig_opts = [
    "-c++",
    "-I/Users/csayres/installTCC/ndarray/8.0.0.0+1/python",
    "-I/Users/csayres/installTCC/coordConv/include",
]


#################### to build coord conv ################
# linkDirs = ["/Users/csayres/installTCC/coordConv/build/lib.macosx-10.7-x86_64-2.7"]
linkDirs = ["/Users/csayres/installTCC/slalib_2013-04-25/install/lib"]
# os.environ["CC"] = "clang++"
module = Extension(
    "python/coordConv/_coordConvLib",
    srcFiles,
    include_dirs=includeList,
    swig_opts=swig_opts,
    library_dirs=linkDirs,
    libraries=['sla'],
    extra_compile_args = extra_compile_args,
    extra_link_args = extra_link_args
)

setup(
    name = "_coordConvLib",
    ext_modules = [module]
)

buildDir = glob.glob("build/lib*")[0]
soFile = glob.glob(buildDir + "/python/coordConv/*so")[0]
base, filename = os.path.split(soFile)
dest = "python/coordConv/%s"%filename
copyfile(soFile, dest)
mode = os.stat(dest).st_mode
mode |= (mode & 0o444) >> 2
os.chmod(dest, mode)

####################################################