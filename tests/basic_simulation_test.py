#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3,
# see "/usr/share/common-licenses/GPL-3".
# Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
# Contributors:
# - Tomasz Makarewicz (makson96@gmail.com)
#
# This tests returs GROMACS version detected by plugin and check if it equal $GROMACS_VER

from __future__ import print_function

import os
import shutil
import sys

import pymol_plugin_dynamics

print("######################")
print("Starting Basic Simulation test")
print("######################")

status, s_params = pymol_plugin_dynamics.init_function(travis_ci=True)

# Download molecule for tests
project_name = "2fjz"
s_params.change_project_name(project_name)
project_dir = pymol_plugin_dynamics.get_project_dirs(project_name)
shutil.copy("/usr/share/pdb-files/2fjz.pdb", project_dir + "2fjz.pdb")

s_params.create_cfg_files()

# Execute dynamics simulation
pymol_plugin_dynamics.dynamics(s_params)

if not os.path.isfile(project_dir + "2fjz_multimodel.pdb"):
	sys.exit(1)

print("######################")
print("Basic Simulation test finished successfully")
print("######################")
