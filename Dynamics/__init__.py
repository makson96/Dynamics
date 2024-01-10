#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3,
# see "/usr/share/common-licenses/GPL-3".

from pymol import plugins
from version import VERSION
import pymol_plugin_dynamics


def __init_plugin__(app=None):
    p_label = f"Dynamics Gromacs {VERSION}"
    plugins.addmenuitemqt(p_label, pymol_plugin_dynamics.init_function)
