#!/usr/bin/env python2
#-*- coding: utf-8 -*-

import sys

sys.path.insert(0,"/usr/lib/pymodules/python2.7/pmg_tk/startup/")
sys.path.insert(0,"/usr/lib/python2.7/dist-packages/pmg_tk/startup/")
sys.path.insert(0,"/usr/local/lib/python2.7/dist-packages/pmg_tk/startup/")

import dynamics_pymol_plugin

dynamics_pymol_plugin.init_function()
