import os
import re
import sys
import platform
import subprocess
from pathlib import Path

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

web_templates = "src/rasputin/web/*"


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON'
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg, f"-DCMAKE_PREFIX_PATH={Path.home() / '.local'};{Path().cwd() / 'lib'}"]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        cmake_build_command = ['cmake', "-S", ".", "-B", "build"] + cmake_args
        _compile_and_check(cmake_build_command, env=env)
        _compile_and_check(["cmake", "--build", "build"])

        print()  # Add an empty line for cleaner output

def _compile_and_check(cmd, *args, **kwargs):
        proc = subprocess.run(cmd, *args, **kwargs, check=True)
        try:
            proc.check_returncode()
        except subprocess.CalledProcessError as err:
            raise RuntimeError(f"{cmd} exited with {proc.returncode}\n{proc.stderr}") from err

setup(
    name='rasputin',
    version='0.1',
    author='Ola Skavhaug',
    author_email='ola@xal.no',
    description='A simple tool for constructing TINs from DEMs',
    long_description='',
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    requires=[
          'Cython',
          'numpy',
          'pyproj',
          'Pillow',
          'h5py',
          'lxml',
          'shapely',
          'descartes',
          'meshio',
    ],
    extras_require={'all': ['pytest']},
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    data_files=[("rasputin/web", ["web/index.js", "web/index.html", "web/data.js", "web/favicon.ico"]),
                ("rasputin/web/js", ["web/js/three.min.js", "web/js/dat.min.js"]),
                ("rasputin/web/js/controls", ["web/js/controls/OrbitControls.js", "web/js/controls/PointerLockControls.js"])],
    # add an extension module named 'python_cpp_example' to the package
    # 'python_cpp_example'
    ext_modules=[CMakeExtension('rasputin/rasputin')],
    # add custom build_ext command
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    entry_points={
        'console_scripts': ['rasputin_web = rasputin.web_visualize:visualize_tin',
                            'rasputin_store = rasputin.application:store_tin',
                            'rasputin_shade = rasputin.application:compute_shades']}
)

