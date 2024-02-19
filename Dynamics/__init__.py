#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3,
# see "/usr/share/common-licenses/GPL-3".

import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.resolve()))

from pymol import plugins

from version import VERSION
import pymol_plugin_dynamics
import GromacsInput
import GromacsOutput
import SimulationParameters
import Vectors
import MdpConfig
import CalculationWindow
def __init_plugin__(app=None):
    p_label = f"Dynamics Gromacs {VERSION}"
    plugins.addmenuitemqt(p_label, pymol_plugin_dynamics.init_function)
