from setuptools import setup, Extension, find_packages
import os
import numpy
import glob
import sys
from shutil import copyfile

coordConvDir = os.getenv("COORDCONV_DIR")
eigenDir = os.getenv("EIGEN_DIR")
boostDir = os.getenv("BOOST_DIR")
ndarrayDir = os.getenv("NDARRAY_DIR")
pythonDir = os.getenv("CONDA_DIR")
condaDir = os.getenv("CONDA_DIR")



srcFiles = [
    coordConvDir + "/python/coordConv/coordConvLib.i",
]

# slalibFiles = glob.glob("/Users/csayres/installTCC/slalib_2013-04-25/src/*.c")

# srcFiles += slalibFiles

convFiles = glob.glob(coordConvDir + "/src/*.cc")
srcFiles += convFiles

includeList = [
    coordConvDir + "/python",
    coordConvDir + "/include",
    numpy.get_include(),
    ndarrayDir + "/python",
    ndarrayDir + "/include",
    coordConvDir + "/include",
    condaDir + "/include/python2.7",
    eigenDir + "/include",
    boostDir,
]

extra_compile_args = ["--std=c++11", "-fPIC", "-v", "-O3"]
extra_link_args = None
if sys.platform == 'darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.9']
    extra_link_args = ["-v", '-mmacosx-version-min=10.9']

swig_opts = [
    "-c++",
    "-I/%s/python"%ndarrayDir,
    "-I%s/include"%coordConvDir
]


#################### to build coord conv ################
# linkDirs = [coordConvDir + "/build/lib.macosx-10.7-x86_64-2.7"]
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