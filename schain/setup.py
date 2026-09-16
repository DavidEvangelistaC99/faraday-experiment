# Copyright (c) 2012-2020 Jicamarca Radio Observatory
# All rights reserved.
#
# Distributed under the terms of the BSD 3-clause license.
"""schainpy is an open source library to read, write and process radar data

Signal Chain is a radar data processing library wich includes modules to read,
and write different files formats, besides modules to process and visualize the
data.
"""

import os
import re
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from schainpy import __version__

DOCLINES = __doc__.split("\n")

gfor = os.popen("gfortran --version").read()
match = re.search(r'\d+',gfor)
gfor_ver = int(match.group())

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(
    name = "schainpy",
    version = __version__,
    description = DOCLINES[0],
    long_description = "\n".join(DOCLINES[2:]),
    url = "https://github.com/JRO-Peru/schainpy",
    author = "Jicamarca Radio Observatory",
    author_email = "jro-developers@jro.igp.gob.pe",
    license="BSD-3-Clause",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
    ],
    packages = {
        'schainpy',
        'schainpy.model',
        'schainpy.model.data',
        'schainpy.model.graphics',
        'schainpy.model.io',
        'schainpy.model.proc',
        'schainpy.model.utils',
        'schainpy.utils',
        'schainpy.gui',
        'schainpy.cli',
        },
    package_data = {'': ['schain.conf.template'],
                    'schainpy.files': ['*.oga']
                    },
    include_package_data = False,
    scripts = ['schainpy/gui/schainGUI'],
    entry_points = {
        'console_scripts': [
            'schain = schainpy.cli.cli:main',
            ],
        },
    cmdclass = {'build_ext': build_ext},
    ext_modules=[
        Extension("schainpy.model.data._noise", ["schainc/_noise.c"]),
        Extension("schainpy.model.data._HS_algorithm", ["schainc/_HS_algorithm.c"]),
        ],
    setup_requires = ["numpy"],
    install_requires = [
        "scipy",
        "h5py",
        "matplotlib",
        "pyzmq",
        "fuzzywuzzy",
        "click",
        ],
)

main_path = os.getcwd()
child_path = '/schainf/Ffiles/bfmodel/'
child_path_r = '/schainf/Ffiles/jlib26feb2001'
#print("main_path: ", main_path)
my_str = '      '+"ppath = "+"\""+main_path+"\"\n"
my_str_2 = '      '+"cpath = "+"\""+child_path+"\""
my_str_2_r = '      '+"cpath = "+"\""+child_path_r+"\""

#'''
with open('./schainf/Ffiles/get_path_1.f', 'r') as path1: data1 = path1.read()
with open('./schainf/Ffiles/get_path_2.f', 'r') as path2: data2 = path2.read()
with open('./schainf/Ffiles/get_path.f', 'w') as final: final.write(data1+my_str+my_str_2+data2)

with open('./schainf/Ffiles/get_path_1_reader.f', 'r') as p1: data1_r = p1.read()
with open('./schainf/Ffiles/get_path_2_reader.f', 'r') as p2: data2_r = p2.read()
with open('./schainf/Ffiles/get_path_reader.f', 'w') as final_r: final_r.write(data1_r+my_str+my_str_2_r+data2_r)
#'''

from numpy.distutils.core import Extension, setup

if gfor_ver >= 10:
    extra_f77 = "-fallow-argument-mismatch"
else:
    extra_f77 = "-g"

setup(name='schainpy',
    ext_modules = [
        Extension("schainpy.model.proc.mkfact_short_2020_2",
            sources=[
                "schainf/Ffiles/mkfact_short_2020_2.pyf",
                "schainf/Ffiles/lmdif1.f",
                "schainf/Ffiles/mkfact.f",
                "schainf/Ffiles/r1mach.f",
                "schainf/Ffiles/bfield2.f",
                "schainf/Ffiles/get_path.f"],
            extra_f77_compile_args=[extra_f77]),
        Extension("schainpy.model.proc.fitacf_guess",
            sources=[
                "schainf/Ffiles/fitacf_guess.pyf",
                "schainf/Ffiles/fitacf_guess.f"],
            extra_f77_compile_args=[extra_f77]),
        Extension("schainpy.model.proc.fitacf_acf2",
            sources = [
                "schainf/Ffiles/fitacf_acf2.pyf",
                "schainf/Ffiles/full_profile_profile.f",
                "schainf/Ffiles/fitacf.f",
                "schainf/Ffiles/get_path_reader.f",
                "schainf/Ffiles/r1mach.f",
                "schainf/Ffiles/lmdif1.f",
                "schainf/Ffiles/lagp.f",
                "schainf/Ffiles/reader.c",
                "schainf/Ffiles/cbesi.f",
                "schainf/Ffiles/i1mach.f",
                "schainf/Ffiles/zeta.f",
                "schainf/Ffiles/qc25f.f",
                "schainf/Ffiles/qwgtf.f",
                "schainf/Ffiles/qcheb.f",
                "schainf/Ffiles/sgtsl.f",
                "schainf/Ffiles/qk15w.f",
                "schainf/Ffiles/complex.c",
                "schainf/Ffiles/cbinu.f",
                "schainf/Ffiles/cseri.f",
                "schainf/Ffiles/cwrsk.f",
                "schainf/Ffiles/crati.f",
                "schainf/Ffiles/casyi.f",
                "schainf/Ffiles/cbuni.f",
                "schainf/Ffiles/cuni2.f",
                "schainf/Ffiles/gamln.f",
                "schainf/Ffiles/cuchk.f",
                "schainf/Ffiles/cbknu.f",
                "schainf/Ffiles/cshch.f",
                "schainf/Ffiles/ckscl.f",
                "schainf/Ffiles/cuoik.f",
                "schainf/Ffiles/cunik.f",
                "schainf/Ffiles/cuni1.f",
                "schainf/Ffiles/cairy.f",
                "schainf/Ffiles/cmlri.f",
                "schainf/Ffiles/cunhj.f",
                "schainf/Ffiles/cacai.f",
                "schainf/Ffiles/csisl.f",
                "schainf/Ffiles/caxpy.f",
                "schainf/Ffiles/cs1s2.f",
                "schainf/Ffiles/scabs1.f",
                "schainf/Ffiles/cdotu.f",
                "schainf/Ffiles/rs.f",
                "schainf/Ffiles/sppfa.f",
                "schainf/Ffiles/sdot.f",
                "schainf/Ffiles/tred2.f",
                "schainf/Ffiles/tql2.f",
                "schainf/Ffiles/sppdi.f",
                "schainf/Ffiles/saxpy.f",
                "schainf/Ffiles/sscal.f",
                "schainf/Ffiles/pythag.f",
                "schainf/Ffiles/tql1.f",
                "schainf/Ffiles/get_path.f",
                "schainf/Ffiles/tred1.f"],
            extra_f77_compile_args=[extra_f77]),
        Extension("schainpy.model.proc.fitacf_fit_short",
            sources = [
                "schainf/Ffiles/fitacf_fit_short.pyf",
                "schainf/Ffiles/cairy.f",
                "schainf/Ffiles/casyi.f",
                "schainf/Ffiles/cbesi.f",
                "schainf/Ffiles/cbinu.f",
                "schainf/Ffiles/cbknu.f",
                "schainf/Ffiles/cbuni.f",
                "schainf/Ffiles/ckscl.f",
                "schainf/Ffiles/crati.f",
                "schainf/Ffiles/cacai.f",
                "schainf/Ffiles/cmlri.f",
                "schainf/Ffiles/cs1s2.f",
                "schainf/Ffiles/cseri.f",
                "schainf/Ffiles/cshch.f",
                "schainf/Ffiles/cuchk.f",
                "schainf/Ffiles/cunhj.f",
                "schainf/Ffiles/cuni1.f",
                "schainf/Ffiles/complex.c",
                "schainf/Ffiles/cuni2.f",
                "schainf/Ffiles/cunik.f",
                "schainf/Ffiles/cuoik.f",
                "schainf/Ffiles/cwrsk.f",
                "schainf/Ffiles/fitacf_fit_short.f",
                "schainf/Ffiles/gamln.f",
                "schainf/Ffiles/i1mach.f",
                "schainf/Ffiles/lmdif1.f",
                "schainf/Ffiles/pythag.f",
                "schainf/Ffiles/qc25f.f",
                "schainf/Ffiles/qcheb.f",
                "schainf/Ffiles/qk15w.f",
                "schainf/Ffiles/qwgtf.f",
                "schainf/Ffiles/r1mach.f",
                "schainf/Ffiles/reader.c",
                "schainf/Ffiles/rs.f",
                "schainf/Ffiles/saxpy.f",
                "schainf/Ffiles/sdot.f",
                "schainf/Ffiles/sgtsl.f",
                "schainf/Ffiles/sppdi.f",
                "schainf/Ffiles/sppfa.f",
                "schainf/Ffiles/sscal.f",
                "schainf/Ffiles/tql1.f",
                "schainf/Ffiles/tql2.f",
                "schainf/Ffiles/tred1.f",
                "schainf/Ffiles/tred2.f",
                "schainf/Ffiles/get_path.f",
                "schainf/Ffiles/get_path_reader.f",
                "schainf/Ffiles/zeta.f"],
            extra_f77_compile_args=[extra_f77]),
        Extension("schainpy.model.proc.full_profile_profile",
            sources = [
                "schainf/Ffiles/full_profile_profile.pyf",
                "schainf/Ffiles/full_profile_profile.f",
                "schainf/Ffiles/get_path_reader.f",
                "schainf/Ffiles/fitacf.f",
                "schainf/Ffiles/r1mach.f",
                "schainf/Ffiles/lmdif1.f",
                "schainf/Ffiles/reader.c",
                "schainf/Ffiles/cbesi.f",
                "schainf/Ffiles/lagp.f",
                "schainf/Ffiles/i1mach.f",
                "schainf/Ffiles/zeta.f",
                "schainf/Ffiles/qc25f.f",
                "schainf/Ffiles/qwgtf.f",
                "schainf/Ffiles/qcheb.f",
                "schainf/Ffiles/sgtsl.f",
                "schainf/Ffiles/qk15w.f",
                "schainf/Ffiles/cbinu.f",
                "schainf/Ffiles/complex.c",
                "schainf/Ffiles/cseri.f",
                "schainf/Ffiles/cwrsk.f",
                "schainf/Ffiles/crati.f",
                "schainf/Ffiles/casyi.f",
                "schainf/Ffiles/cbuni.f",
                "schainf/Ffiles/cuni2.f",
                "schainf/Ffiles/gamln.f",
                "schainf/Ffiles/cuchk.f",
                "schainf/Ffiles/cbknu.f",
                "schainf/Ffiles/cshch.f",
                "schainf/Ffiles/ckscl.f",
                "schainf/Ffiles/cuoik.f",
                "schainf/Ffiles/cunik.f",
                "schainf/Ffiles/cuni1.f",
                "schainf/Ffiles/cairy.f",
                "schainf/Ffiles/cmlri.f",
                "schainf/Ffiles/cunhj.f",
                "schainf/Ffiles/cacai.f",
                "schainf/Ffiles/csisl.f",
                "schainf/Ffiles/caxpy.f",
                "schainf/Ffiles/cs1s2.f",
                "schainf/Ffiles/scabs1.f",
                "schainf/Ffiles/cdotu.f",
                "schainf/Ffiles/rs.f",
                "schainf/Ffiles/sppfa.f",
                "schainf/Ffiles/sdot.f",
                "schainf/Ffiles/tred2.f",
                "schainf/Ffiles/tql2.f",
                "schainf/Ffiles/sppdi.f",
                "schainf/Ffiles/saxpy.f",
                "schainf/Ffiles/sscal.f",
                "schainf/Ffiles/pythag.f",
                "schainf/Ffiles/tql1.f",
                "schainf/Ffiles/get_path.f",
                "schainf/Ffiles/tred1.f"],
            extra_f77_compile_args=[extra_f77])
                ]
                )

