#!/bin/bash

# This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
# Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
# Contributors:
# - Tomasz Makarewicz (makson96@gmail.com)
#

# Note: this script should be run from main Dynamics PyMOL Plugin Directory

export PYTHONPATH=$PWD:$PYTHONPATH

if ! python3 tests/unit_tests.py ; then
exit 1
fi

if ! python3 tests/init_test.py ; then
exit 1
fi

if ! python3 tests/gromacs_ver_test.py ; then
exit 1
fi

if ! python3 tests/basic_simulation_test.py ; then
exit 1
fi
