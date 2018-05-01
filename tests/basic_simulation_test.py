#!/usr/bin/env python2
#-*- coding: utf-8 -*-

##This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
##Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
##Contributors:
##- Tomasz Makarewicz (makson96@gmail.com)
##
#This tests returs GROMACS version detected by plugin and check if it equal $GROMACS_VER

print "######################"
print "Starting Basic Simulation test"
print "######################"

import os, sys
one_up = os.path.os.path.abspath(__file__)
one_up = one_up.split("tests")[0]
sys.path.insert(0, one_up)

import pymol_plugin_dynamics
pymol_plugin_dynamics.init_function(travisCI=True)

#Download molecule for tests
os.makedirs(pymol_plugin_dynamics.project_dir)
import urllib
urllib.urlretrieve ("https://files.rcsb.org/download/1RO3.pdb", pymol_plugin_dynamics.project_dir + "1ro3.pdb")

status = ["ok", "ok"]
stop = 0
project_name = "1ro3"

pymol_plugin_dynamics.dynamics(pymol_plugin_dynamics.project_dir + "1ro3_multimodel.pdb") == False:
	sys.exit(1)

if os.path.isfile()

print "######################"
print "Basic Simulation test finished successfully"
print "######################"
