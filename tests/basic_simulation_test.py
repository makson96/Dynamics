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

pymol_plugin_dynamics.init_function(travis_ci=True)

# Download molecule for tests
pymol_plugin_dynamics.project_name = "2fjz"
pymol_plugin_dynamics.set_gromacsProject_dir()
os.makedirs(pymol_plugin_dynamics.project_dir)
shutil.copy("/usr/share/pdb-files/2fjz.pdb", pymol_plugin_dynamics.project_dir + "2fjz.pdb")

pymol_plugin_dynamics.status = ["ok", "ok"]
pymol_plugin_dynamics.stop = 0
pymol_plugin_dynamics.prody_true = 0
pymol_plugin_dynamics.create_config_files()

# Execute dynamics simulation
pymol_plugin_dynamics.dynamics()

if not os.path.isfile(pymol_plugin_dynamics.project_dir + "2fjz_multimodel.pdb"):
	sys.exit(1)

print("######################")
print("Basic Simulation test finished successfully")
print("######################")
