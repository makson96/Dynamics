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

import os, sys, shutil
one_up = os.path.os.path.abspath(__file__)
one_up = one_up.split("tests")[0]
sys.path.insert(0, one_up)

import pymol_plugin_dynamics
pymol_plugin_dynamics.init_function(travisCI=True)

#Download molecule for tests
pymol_plugin_dynamics.project_name = "2fjz"
pymol_plugin_dynamics.setGromacsProjectDir()
os.makedirs(pymol_plugin_dynamics.project_dir)
shutil.copy("/usr/share/pdb-files/2fjz.pdb", pymol_plugin_dynamics.project_dir + "2fjz.pdb")

pymol_plugin_dynamics.status = ["ok", "ok"]
pymol_plugin_dynamics.stop = 0
pymol_plugin_dynamics.prody_true = 0
pymol_plugin_dynamics.create_config_files()

#Execute dynamics simulation
pymol_plugin_dynamics.dynamics()

if os.path.isfile(pymol_plugin_dynamics.project_dir + "2fjz_multimodel.pdb") == False:
	sys.exit(1)

print "######################"
print "Basic Simulation test finished successfully"
print "######################"
