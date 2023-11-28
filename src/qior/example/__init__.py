# Copyright 2023 and later, Andres Agusti Casado
# This file is part of the python package qior.
# qior is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
# qior is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
# You should have received a copy of the GNU General Public License along 
# with qior. If not, see <https://www.gnu.org/licenses/>. 
"""
This is a script that copies a bunch of scripts to the directory qior is
being executed from. Those scripts are examples showcasing qior's features.
"""

import os, shutil

def main():
    deploy_examples()

def deploy_examples():
    print("copying examples to %s" % user_directory())
    print("run them with 'python exampleX.py' with 'X' an integer character, and in an environment with 'qior' installed")
    for filename in get_example_filenames():
        copy_to(filename, user_directory())

def user_directory():
    return os.path.realpath(os.path.curdir)

def get_example_filenames():
    directory = os.path.dirname(os.path.realpath(__file__))
    return [os.path.join(directory, fn) for fn in os.listdir(directory) if "example" in fn]


def copy_to(filename, directory):
    shutil.copy(filename, directory)
