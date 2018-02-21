#!/usr/bin/env python2
#-*- coding: utf-8 -*-

##This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
##Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
##Contributors:
##- Tomasz Makarewicz (makson96@gmail.com)
##

import os, sys
one_up = os.path.os.path.abspath(__file__)
one_up = one_up.split("tests")[0]
sys.path.insert(0, one_up)

print "This is simple initialization check script"
import pymol_plugin_dynamics
pymol_plugin_dynamics.init_function(travisCI=True)
