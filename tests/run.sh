##This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
##Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
##Contributors:
##- Tomasz Makarewicz (makson96@gmail.com)
##

#Note: this script should be run from main Dynamics PyMOL Plugin Directory

if ! python tests/init_test.py ; then
exit 1
fi

if ! python tests/gromacs_ver_test.py ; then
exit 1
fi
