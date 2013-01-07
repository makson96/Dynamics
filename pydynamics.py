#!/usr/bin/env python2
#-*- coding: utf-8 -*-

import sys, os, shutil

current_dir = os.getcwd()
sys.path.insert(0,"/usr/lib/pymodules/python2.7/pmg_tk/startup/")
sys.path.insert(0,"/usr/lib/python2.7/dist-packages/pmg_tk/startup/")
sys.path.insert(0,"/usr/local/lib/python2.7/dist-packages/pmg_tk/startup/")

import dynamics_pymol_plugin
		
if os.path.isfile(sys.argv[1]) and len(sys.argv) == 3:
	project_name, project_dir = dynamics_pymol_plugin.dynamics_cmd(sys.argv[1], [sys.argv[1]])
	print project_dir + project_name + "_multimodel.pdb"
	print sys.argv[2]
	os.chdir(current_dir)
	shutil.copy(project_dir + project_name + "_multimodel.pdb", sys.argv[2])
else:
	help_name, clean_name = dynamics_pymol_plugin.init_function(1)
	if help_name.count(sys.argv[1]) == 1 or clean_name.count(sys.argv[1]) == 1:
		dynamics_pymol_plugin.dynamics(sys.argv[1])
	else:
		print "Wrong usage of PyDynamics. Try 'pydynamics -h' to get some help informations"
