#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3, see
# "/usr/share/common-licenses/GPL-3".
# Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
# Contributors:
# - Tomasz Makarewicz (makson96@gmail.com)
#
# This tests returs GROMACS version detected by plugin and check if it equal $GROMACS_VER

from __future__ import print_function

import os
import sys

import pymol_plugin_dynamics

print("######################")
print("Starting GROMACS version test")
print("######################")

pymol_plugin_dynamics.init_function(travis_ci=True)
print("Detected GROMACS version is: " + pymol_plugin_dynamics.gmx_version)
if os.environ['GROMACS_VER'] != pymol_plugin_dynamics.gmx_version:
	log = open(os.environ['HOME'] + "/.dynamics/test_gromacs.txt", "r").read()
	print(log)
	sys.exit(1)
	
print("######################")
print("GROMACS version test finished successfully")
print("######################")
