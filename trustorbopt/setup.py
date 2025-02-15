import os
import sys
import subprocess
from setuptools import setup, find_packages
from setuptools.command.build import build


ext = "dylib" if sys.platform == "darwin" else "so"


class CMakeBuild(build):
    def run(self):
        # create build directory
        build_dir = os.path.abspath("build")
        os.makedirs(build_dir, exist_ok=True)

        # run CMake configuration step
        subprocess.check_call(["cmake", ".."], cwd=build_dir)

        # run CMake build step
        subprocess.check_call(["cmake", "--build", "."], cwd=build_dir)

        # copy the built library into the package directory
        lib_file = "libtrustorbopt." + ext
        src_lib_path = os.path.join(build_dir, lib_file)
        dest_lib_path = os.path.join(os.path.abspath("pytrustorbopt"), lib_file)
        subprocess.check_call(["cp", src_lib_path, dest_lib_path])

        # run steps in parent class
        super().run()


setup(
    name="pytrustorbopt",
    version="0.1",
    author="Jonas Greiner",
    author_email="jongr@kemi.dtu.dk",
    description="Python interface to trustorbopt Fortran library",
    packages=find_packages(),
    package_data={"pytrustorbopt": ["libtrustorbopt." + ext]},
    include_package_data=True,
    cmdclass={"build": CMakeBuild},
    install_requires=["numpy"],
)
