#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3,
# see "/usr/share/common-licenses/GPL-3".
# Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
# Contributors:
# - Tomasz Makarewicz (makson96@gmail.com)
#
# This is simple initialization check script

from __future__ import print_function

import os
import sys

import pymol_plugin_dynamics

print("######################")
print("Starting Initialization test")
print("######################")

pymol_plugin_dynamics.init_function(travis_ci=True)

print("######################")
print("Initialization finished successfully")
print("######################")
