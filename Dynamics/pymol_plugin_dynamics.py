#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This software (including its Debian packaging) is available to you under the terms of the GPL-3,
# see "/usr/share/common-licenses/GPL-3".
# Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
# Contributors:
# - Tomasz Makarewicz (makson96@gmail.com)
# - Ajit B. Datta (ajit@jcbose.ac.in)
# - Sara Boch Kminikowska
# - Manish Sud (msud@san.rr.com; URL: www.MayaChemTools.org)
# - Thomas Holder
#

from __future__ import print_function

# --Import libraries--
# Import nativ python libraries
import os
import pickle
import shutil
import subprocess
import sys
import time
import tarfile
# This is actually needed if Tk will be removed.
import re

# Import libraries for tk graphic interface
import _thread as thread
import queue as Queue
from tkinter import *
from tkinter import messagebox as tkMessageBox
from tkinter import filedialog as tkFileDialog
from tkinter.ttk import Progressbar, Scrollbar
# Import libraries from PyMOL specific work.
from pymol import cmd, cgo, parsing, plugins, CmdException

# TODO: It seams that stored is removed from PyMOL API. We need to handle it correctly
try:
    from pymol import stored
except ImportError:
    stored = False

# Check for ProDy
try:
    import prody
except ModuleNotFoundError:
    prody = False

# Plugin Version
from version import VERSION
plugin_ver = f" {VERSION}"

EM_INIT_CONFIG = """define = -DFLEX_SPC
constraints = none
integrator = steep
nsteps = 10000
nstlist = 10
ns_type = simple
rlist = 1.5
rcoulomb = 1.5
rvdw = 1.5
emtol = 1000.0
emstep = 0.01
implicit-solvent = no
;gb-algorithm = Still
;pbc = no
;rgbradii = 0
cutoff-scheme = Verlet
coulombtype = PME"""

PR_INIT_CONFIG = """define = -DPOSRES
constraints = all-bonds
integrator = md-vv
dt = 0.002
nsteps = 5000
nstcomm = 1
nstxout = 100
nstvout = 100
nstfout = 0
nstlog = 10
nstenergy = 10
nstlist = 10
ns_type = simple
rlist = 1.5
rcoulomb = 1.5
rvdw = 1.5
Tcoupl = v-rescale
tau_t = 0.1 0.1
tc-grps = protein Non-Protein
ref_t = 298 298
Pcoupl = no
tau_p = 0.5
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = yes
gen_temp = 298.0
gen_seed = 173529
cutoff-scheme = Verlet
coulombtype = PME"""

MD_INIT_CONFIG = """;define = -DPOSRES
integrator = md-vv
dt = 0.002
nsteps = 5000
nstcomm = 1
nstxout = 50
nstvout = 50
nstfout = 0
nstlist = 10
ns_type = simple
rlist = 1.5
rcoulomb = 1.5
rvdw = 1.5
Tcoupl = v-rescale
tau_t = 0.1 0.1
tc-grps = protein Non-Protein
ref_t = 298 298
Pcoupl = no
tau_p = 0.5
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = yes
gen_temp = 298.0
gen_seed = 173529
constraints = all-bonds
constraint-algorithm = Lincs
continuation = no
shake-tol = 0.0001
lincs-order = 4
lincs-warnangle = 30
morse = no
implicit-solvent = no
;gb-algorithm = Still
;pbc = no
;rgbradii = 0
;comm_mode = ANGULAR
cutoff-scheme = Verlet
coulombtype = PME"""


# This function will initialize all plugin stufs
def init_function(travis_ci=False, gui_library="qt", parent=False):
    # Fallback to tk, till qt is ready
    gui_library = "tk"
    status = ["ok", ""]
    # Make sure HOME environment variable is defined before setting up directories...
    home_dir = os.path.expanduser('~')
    if home_dir:
        os.chdir(home_dir)
    else:
        print("HOME environment variable not defined")
        status = ["fail", "HOME environment variable not defined. Please set its value and try again."]
    dynamics_dir = get_dynamics_dir()
    project_dir = get_project_dirs()
    # Clean up any temporary project directory...
    if os.path.isdir(project_dir):
        shutil.rmtree(project_dir)
    # Create temporary project directory along with any subdirectories...
    if not os.path.isdir(project_dir):
        os.makedirs(project_dir)

    print("Searching for GROMACS installation")
    os.chdir(dynamics_dir)
    gmx_exe, gmx_version, gmx_build_arch, gmx_on_cygwin = get_gromacs_exe_info()
    os.chdir(home_dir)

    supported_gmx_versions = ["2016", "2018"]
    if not len(gmx_exe):
        print("GROMACS 2016 or newer not detected.")
        status = ["fail",
                  "GROMACS not detected. Please install and setup GROMACS 2016 or newer correctly for your platform."
                  " Check '~/.dynamics/test_gromacs.txt' for more details. Don't forget to add GROMACS bin directory"
                  " to your PATH"]
    elif gmx_version[0:4] not in supported_gmx_versions:
        print("Warning. Unsupported GROMACS Version")
    if status[0] == "ok":
        simulation_parameters = SimulationParameters()
    else:
        simulation_parameters = False
    if not travis_ci:
        create_gui(gui_library, status, simulation_parameters, parent)
    return status, simulation_parameters


class SimulationParameters:
    gmx_output = ""
    gmx_input = ""
    vectors_prody = False
    stop = False
    project_name = "nothing"
    progress = ""
    em_file = ""
    pr_file = ""
    md_file = ""

    def __init__(self):
        self.gmx_output = GromacsOutput()
        self.gmx_input = GromacsInput()
        print("Found GROMACS VERSION {}".format(self.gmx_output.version))
        if prody:
            self.vectors_prody = Vectors()
            print("ProDy correctly imported")
        self.progress = ProgressStatus()

    def create_cfg_files(self):
        self.em_file, self.pr_file, self.md_file = create_config_files(self.project_name)

    def change_stop_value(self, value):
        if value:
            self.stop = True
        else:
            self.stop = False

    def change_project_name(self, name):
        self.project_name = name
        project_dir = get_project_dirs(self.project_name)
        if not os.path.isdir(project_dir):
            os.makedirs(project_dir)


# This class is responsible for interface to GROMACS. It will read all important data from GROMACS tools.
class GromacsOutput:
    version = "GROMACS not found"
    command = ""
    force_list = []
    water_list = []
    group_list = []
    restraints = []

    def __init__(self):
        # Remove garbage
        dynamics_dir = get_dynamics_dir()
        garbage_files = next(os.walk(dynamics_dir))[2]
        for garbage in garbage_files:
            if garbage[0] == "#":
                os.remove(dynamics_dir + garbage)

        gmx_exe, gmx_version, gmx_build_arch, gmx_on_cygwin = get_gromacs_exe_info()
        self.version = gmx_version
        self.command = gmx_exe

        # Track current directiry and switch to dynamics_dir before invoking gmx...
        current_dir = os.getcwd()
        os.chdir(dynamics_dir)

        self.init2()

        # Switch back to current directory...
        os.chdir(current_dir)

    def init2(self):
        print("Reading available force fields and water models")

        fo = open("test_gromacs.pdb", "wb")
        fo.write(b"ATOM      1  N   LYS     1      24.966  -0.646  22.314  1.00 32.74      1SRN  99\n")
        fo.close()

        gmx_stdin_file_path = "gromacs_stdin.txt"
        fo = open(gmx_stdin_file_path, "w")
        fo.write("1\n")
        fo.write("1")
        fo.close()

        gmx_stdout_file_path = "test_gromacs.txt"

        cmd = "{} pdb2gmx -f test_gromacs.pdb -o test_gromacs.gro -p test_gromacs.top".format(self.command)
        execute_subprocess(cmd, gmx_stdin_file_path, gmx_stdout_file_path)
        lista_gromacs = read_text_lines(gmx_stdout_file_path)

        # Reading available force fields
        force_start_line = 0
        while lista_gromacs[force_start_line] != "Select the Force Field:\n":
            force_start_line = force_start_line + 1
        force_start_line = force_start_line + 2
        force_end_line = force_start_line
        while lista_gromacs[force_end_line] != "\n":
            force_end_line = force_end_line + 1
        force_list = lista_gromacs[force_start_line:force_end_line]
        force_list2 = []
        number = 1
        for force in force_list:
            force_list2.append([number, force[:-1]])
            number = number + 1

        self.force_list = force_list2

        # Reading available water models
        self.water_list = get_water_models_info(lista_gromacs)

        print("Reading available groups")

        gmx_stdin_file_path = "gromacs_stdin.txt"
        fo = open(gmx_stdin_file_path, "w")
        fo.write("1")
        fo.close()

        gmx_stdout_file_path = "test_gromacs.txt"

        cmd = "{} trjconv -f  test_gromacs.pdb -s test_gromacs.pdb -o test_gromacs2.pdb".format(self.command)
        execute_subprocess(cmd, gmx_stdin_file_path, gmx_stdout_file_path)

        group_test_list = read_text_lines(gmx_stdout_file_path)

        # Reading available groups
        group_start_line = 0
        while group_test_list[group_start_line] != "Will write pdb: Protein data bank file\n":
            group_start_line = group_start_line + 1
        group_start_line = group_start_line + 1
        group_end_line = group_start_line + 1
        while group_test_list[group_end_line][0:14] != "Select a group":
            group_end_line = group_end_line + 1
        group_list = group_test_list[group_start_line:group_end_line]
        group_list2 = []
        number = 0
        for group in group_list:
            group1 = group.split(' has')
            group2 = group1[0].split('Group     ')
            if len(group2) == 2:
                group_list2.append([number, group2[1]])
                number = number + 1

        self.group_list = group_list2

    # This function will update water list if force field is changed.
    def water_update(self, force_number):
        # Track current directiry and switch to dynamics_dir before invoking gmx...
        current_dir = os.getcwd()
        os.chdir(get_dynamics_dir())

        print("Updating available water models")
        gmx_stdin_file_path = "gromacs_stdin.txt"
        fo = open(gmx_stdin_file_path, "w")
        fo.write("%d\n" % force_number)
        fo.write("1")
        fo.close()

        gmx_stdout_file_path = "test_gromacs.txt"

        cmd = "{} pdb2gmx -f  test_gromacs.pdb -o test_gromacs.gro -p test_gromacs.top".format(self.command)
        execute_subprocess(cmd, gmx_stdin_file_path, gmx_stdout_file_path)

        lista_gromacs = read_text_lines(gmx_stdout_file_path)
        self.water_list = get_water_models_info(lista_gromacs)

        # Switch back to current directory...
        os.chdir(current_dir)

        # save_options()
        return self.water_list

    # This function will read atoms group for restraints for current molecule.
    def restraints_index(self, project_name):
        self.restraints = []
        current_dir = os.getcwd()
        os.chdir(get_project_dirs(project_name))

        fo = open("gromacs_stdin.txt", "w")
        fo.write("q")
        fo.close()

        cmd = "{} make_ndx -f {}.pdb -o index.ndx".format(self.command, project_name)
        execute_subprocess(cmd, "gromacs_stdin.txt", "restraints.log")

        index_list = read_text_lines("restraints.log")

        index_position = 0
        atoms = ""
        for line in index_list:
            if line[0] == "[":
                self.restraints.append([])
                self.restraints[index_position].append(line)
                if index_position != 0:
                    self.restraints[index_position - 1].append(atoms)
                index_position = index_position + 1
                atoms = ""
            else:
                atoms = atoms + line
        self.restraints[index_position - 1].append(atoms)
        os.chdir(current_dir)


# This class is responsible for performing molecular dynamics simulation with GROMACS tools.
class GromacsInput:
    force = 1
    water = 1
    group = 1
    box_type = "triclinic"
    explicit = 1
    # variable to choose heavy hydrogen
    hydro = "noheavyh"
    box_distance = "0.8"
    box_density = "1000"
    restraints_nr = 1
    # four variables salt, positive, negative and neutral
    neutrality = "neutral"
    salt_conc = "0.15"
    positive_ion = "NA"
    negative_ion = "CL"
    command_distinction = "\n!************************!\n"

    # This function will change given variabless stored by the class (needed for lambda statements)
    def update(self, gmx_options):
        for key, value in gmx_options.items():
            if key == "force":
                self.force = value
            elif key == "water":
                self.water = value
            elif key == "group":
                self.group = value
            elif key == "box_type":
                self.box_type = value
            elif key == "hydro":
                self.hydro = value
            elif key == "box_distance":
                self.box_distance = value
            elif key == "box_density":
                self.box_density = value
            elif key == "restraints_nr":
                self.restraints_nr = value
            elif key == "neutrality":
                self.neutrality = value
            elif key == "salt_conc":
                self.salt_conc = value
            elif key == "positive_ion":
                self.positive_ion = value
            elif key == "negative_ion":
                self.negative_ion = value
            elif key == "explicit":
                self.explicit = value
        #        save_options()
        print("gromacs updated")

    # This function will create initial topology and trajectory using pdb file and choosen force field
    def pdb2top(self, s_params):
        status = ["ok", "Calculating topology using Force fields"]
        status_update(status)
        hh = "-" + self.hydro
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}.gro".format(project_name))
            os.remove("{}.top".format(project_name))
        except FileNotFoundError:
            pass

        fo = open("gromacs_stdin.txt", "w")
        fo.write("%s\n" % str(self.force))
        fo.write("%s" % str(self.water))
        fo.close()

        command = "{0} pdb2gmx -f {1}.pdb -o {1}.gro -p {1}.top {2}".format(gmx_cmd, project_name, hh)
        execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        if os.path.isfile("{}.gro".format(project_name)):
            status = ["ok", ""]
        else:
            status = ["fail", "Warning. Trying to ignore unnecessary hydrogen atoms."]

            command = "{0} pdb2gmx -ignh -f {1}.pdb -o {1}.gro -p {1}.top {2}".format(gmx_cmd, project_name, hh)
            execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')
        status_update(status)

        stop = s_params.stop
        if os.path.isfile("{}.gro".format(project_name)) and not stop:
            status = ["ok", "Calculated topology using Force fields"]
        else:
            status = ["fail", "Force field unable to create topology file"]
        return status

    # This is alternative function to create initial topology and triectory using pdb file
    @staticmethod
    def x2top(s_params):
        status = ["ok", "Calculating topology using Force fields"]
        status_update(status)
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}.gro".format(project_name))
            os.remove("{}.top".format(project_name))
        except FileNotFoundError:
            pass

        command = "{0} x2top -f {1}.pdb -o {1}.top".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}.top".format(project_name)) and not stop:
            status = ["ok", "Calculating structure using trjconv."]
        else:
            status = ["fail", "Unable to create topology file."]
        status_update(status)

        if status[0] == "ok":
            fo = open("gromacs_stdin.txt", "w")
            fo.write("0")
            fo.close()

            command = "{0} trjconv -f {1}.pdb -s {1}.pdb -o {1}.gro".format(gmx_cmd, project_name)
            execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

            stop = s_params.stop
            if os.path.isfile("{}.gro".format(project_name)) and not stop:
                status = ["ok", "Calculated structure using trjconv."]
            else:
                status = ["fail", "Unable to create structure file."]
        return status

    # This function will create and add waterbox.
    def waterbox(self, s_params):
        status = ["ok", "Generating waterbox"]
        box_type = "-bt {}".format(self.box_type)
        distance = "-d {}".format(self.box_distance)
        density = "-density {}".format(self.box_density)
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}1.gro".format(project_name))
            os.remove("{}_solv.gro".format(project_name))
        except FileNotFoundError:
            pass

        status_update(status)
        command = "{0} editconf -f {1}.gro -o {1}1.gro -c {2} {3} {4}".format(gmx_cmd, project_name, box_type, distance,
                                                                              density)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        water_name = s_params.gmx_output.water_list[self.water - 1][1][4:8].lower()
        print(water_name)
        if water_name == "tip4":
            water_gro = "tip4p.gro"
        elif water_name == "tip5":
            water_gro = "tip5p.gro"
        else:
            water_gro = "spc216.gro"

        command = "{0} solvate -cp {1}1.gro -cs {2} -o {1}_solv.gro -p {1}.top".format(gmx_cmd, project_name, water_gro)

        status = ["ok", "Adding Water Box"]
        status_update(status)

        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}1.gro".format(project_name)) and not stop:
            status = ["ok", "Water Box Added"]
        else:
            status = ["fail", "Unable to add water box"]
        return status

    # This function will add ions/salts to the protein in waterbox
    def saltadd(self, s_params):
        status = ["ok", "Preparing to add ions or salt"]
        salt = "-conc {}".format(self.salt_conc)
        positive = "-pname {}".format(self.positive_ion)
        negative = "-nname {}".format(self.negative_ion)
        neu = "-{}".format(self.neutrality)
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove(project_name + "_b4em.gro")
            os.remove(project_name + "_ions.tpr")
        except FileNotFoundError:
            pass

        command = "{0} grompp -f em -c {1}_solv.gro -o {1}_ions.tpr -p {1}.top -maxwarn 1".format(gmx_cmd, project_name)
        status_update(status)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        fo = open("gromacs_stdin.txt", "w")
        fo.write("13")
        fo.close()

        status = ["ok", "Adding salts and ions"]
        status_update(status)

        command = "{0} genion -s {1}_ions.tpr -o {1}_b4em.gro {2} {3} {4} {5} -p {1}.top".format(gmx_cmd,
                                                                                                 project_name, positive,
                                                                                                 negative, salt, neu)
        execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_b4em.gro".format(project_name)) and not stop:
            status = ["ok", "Ions added successfully"]
        elif stop == 0:
            status = ["ok", "Find out what's wrong!"]
        else:
            status = ["failed", "Unable to add ions"]
        return status

    # This function will perform energy minimization
    @staticmethod
    def em(s_params):
        status = ["ok", "Energy Minimization"]
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}_em.tpr".format(project_name))
            os.remove("{}_em.trr".format(project_name))
            os.remove("{}_b4pr.gro".format(project_name))
        except FileNotFoundError:
            pass

        # Check if waterbox was added and adjust accordingly.
        if not os.path.isfile("{}_b4em.gro".format(project_name)):
            if os.path.isfile("{}_solv.gro".format(project_name)):
                shutil.copy("{}_solv.gro".format(project_name), "{}_b4em.gro".format(project_name))
            elif os.path.isfile(project_name + "{}.gro".format(project_name)):
                shutil.copy("{}.gro".format(project_name), "{}_b4em.gro".format(project_name))

        status_update(status)
        command = "{0} grompp -f em -c {1}_b4em -p {1} -o {1}_em".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_em -o {1}_em -c {1}_b4pr -v".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_em.tpr".format(project_name)) and os.path.isfile("{}_b4pr.gro".format(project_name)) and \
                not stop:
            status = ["ok", "Energy Minimized"]
        else:
            status = ["fail", "Unable to perform Energy Minimization"]
        return status

    # This function will perform position restrained MD
    @staticmethod
    def pr(s_params):
        status = ["ok", "Position Restrained MD"]
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}_pr.tpr".format(project_name))
            os.remove("{}_pr.trr".format(project_name))
            os.remove("{}_b4md.gro".format(project_name))
        except FileNotFoundError:
            pass

        status_update(status)
        command = "{0} grompp -f pr -c {1}_b4pr -r {1}_b4pr -p {1} -o {1}_pr".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_pr -o {1}_pr -c {1}_b4md -v".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_pr.tpr".format(project_name)) and not stop:
            status = ["ok", "Position Restrained MD finished"]
        else:
            status = ["fail", "Unable to perform Position Restrained"]
        return status

    # This function will create posre.itp file for molecular dynamics simulation with choosen atoms if
    # restraints were selected
    @staticmethod
    def restraints(s_params):
        status = ["ok", "Adding Restraints"]
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("posre_2.itp")
        except FileNotFoundError:
            pass

        fo = open("gromacs_stdin.txt", "w")
        fo.write("0")
        fo.close()

        status_update(status)
        command = "{0} genrestr -f {1}.pdb -o posre_2.itp -n index_dynamics.ndx".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("posre_2.itp") and not stop:
            status = ["ok", "Added Restraints"]
            if os.path.isfile("posre.itp"):
                os.remove("posre.itp")
            shutil.copy("posre_2.itp", "posre.itp")
        else:
            status = ["fail", "Unable to create restraints file"]
        return status

    # This function will perform position final molecular dynamics simulation
    @staticmethod
    def md(s_params):
        status = ["ok", "Molecular Dynamics Simulation"]
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}_md.tpr".format(project_name))
            os.remove("{}_md.trr".format(project_name))
        except FileNotFoundError:
            pass

        # Check if em and/or pr was done and adjust accordingly.
        if not os.path.isfile("{}_b4md.gro".format(project_name)):
            if not os.path.isfile("{}_b4pr.gro".format(project_name)):
                # No em and pr
                shutil.copy("{}_b4em.gro".format(project_name), "{}_b4md.gro".format(project_name))
            else:
                # No pr
                shutil.copy("{}_b4pr.gro".format(project_name), "{}_b4md.gro".format(project_name))

        status_update(status)
        command = "{0} grompp -f md -c {1}_b4md  -p {1} -o {1}_md".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_md -o {1}_md -c {1}_after_md -v".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_md.tpr".format(project_name)) and not stop:
            status = ["ok", "Molecular Dynamics Simulation finished"]
        else:
            status = ["fail", "Unable to perform Molecular Dynamics Simulation"]
        return status

    # This function will convert final results to multimodel pdb file
    def trjconv(self, s_params):
        status = ["ok", "Creating Multimodel PDB"]
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}_multimodel.pdb".format(project_name))
        except FileNotFoundError:
            pass
        if os.path.isfile("{}_multimodel.pdb".format(project_name)):
            os.remove("{}_multimodel.pdb".format(project_name))

        fo = open("gromacs_stdin.txt", "w")
        fo.write("%s" % str(self.group))
        fo.close()

        status_update(status)
        command = "{0} trjconv -f {1}_md.trr -s {1}_md.tpr -o {1}_multimodel.pdb".format(gmx_cmd, project_name)
        execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_multimodel.pdb".format(project_name)) and not stop:
            status = ["ok", "Finished!"]
        else:
            status = ["fail", "Unable to generate multimodel PDB file"]
        return status


# This class will handle PCA by ProDy python library and show vectors from NMD file.
class Vectors:
    nmd_name = []
    nmd_atomnames = []
    nmd_resnames = []
    nmd_resids = []
    nmd_bfactors = []
    nmd_coordinates = []
    nmd_mode = []
    nmd_scale_mode = []
    color = "grey"
    scale = 1.0
    mode_nr = 0

    calculation_type = 0
    contact_map = 0
    block_contact_map = 0

    enm = 0

    # Change Multimodel PDB file into NMD vector file
    def prody(self, project_name):
        # Silence ProDy and create logs
        prody.confProDy(verbosity='none')
        prody.startLogfile("log_prody.log")
        # Prepare ensemble
        model = prody.parsePDB(project_name + "_multimodel.pdb", subset='calpha')
        ensemble = prody.Ensemble(project_name + ' ensemble')
        ensemble.setCoords(model.getCoords())
        ensemble.addCoordset(model.getCoordsets())
        ensemble.iterpose()
        # ANM calculations
        if self.calculation_type == 0:
            anm = prody.ANM(project_name)
            anm.buildHessian(ensemble)
            anm.calcModes()
            write_nmd = anm
            self.enm = anm
        # PCA calculations
        elif self.calculation_type == 1:
            pca = prody.PCA(project_name)
            pca.buildCovariance(ensemble)
            pca.calcModes()
            write_nmd = pca
        # GNM calculations
        elif self.calculation_type == 2:
            gnm = prody.GNM(project_name)
            gnm.buildKirchhoff(ensemble)
            gnm.calcModes()
            write_nmd = gnm
            self.enm = gnm
        # Write NMD file
        prody.writeNMD(project_name + '.nmd', write_nmd[:3], model)
        prody.closeLogfile("log_prody.log")

    # Read NMD file
    def nmd_format(self, project_name):
        file_nmd = open('{}.nmd'.format(project_name), "r")
        list_nmd = file_nmd.readlines()

        self.nmd_mode = []
        self.nmd_scale_mode = []
        for line in list_nmd:
            split_line = line.split()
            if split_line[0] == "name":
                self.nmd_name = split_line
                self.nmd_name.pop(0)
            elif split_line[0] == "atomnames":
                self.nmd_atomnames = split_line
                self.nmd_atomnames.pop(0)
            elif split_line[0] == "resnames":
                self.nmd_resnames = split_line
                self.nmd_resnames.pop(0)
            elif split_line[0] == "resids":
                self.nmd_resids = split_line
                self.nmd_resids.pop(0)
            elif split_line[0] == "bfactors":
                self.nmd_bfactors = split_line
                self.nmd_bfactors.pop(0)
            elif split_line[0] == "coordinates":
                self.nmd_coordinates = split_line
                self.nmd_coordinates.pop(0)
            elif split_line[0] == "mode":
                pre_mode = split_line
                self.nmd_mode.append(pre_mode[3:])
                self.nmd_scale_mode.append(pre_mode[2])

    # Show contact map on PyMOL screen
    def show_contact_map(self, sensitivity, project_name):
        contact_matrix = self.enm.getKirchhoff()
        print(contact_matrix)
        c_alpha_nr = 0
        for c_alpha_list in contact_matrix:
            c_alpha_nr = c_alpha_nr + 1
            c_alpha_target_nr = 0
            for c_alpha_1 in c_alpha_list:
                c_alpha_target_nr = c_alpha_target_nr + 1
                if c_alpha_nr != c_alpha_target_nr and float(c_alpha_1) < float(sensitivity):
                    cmd.select("sele1",
                               "n. ca and {}_multimodel and i. {}".format(project_name, str(c_alpha_nr)))  # PyMOL API
                    cmd.select("sele2", "n. ca and {}_multimodel and i. {}".format(project_name,
                                                                                   str(c_alpha_target_nr)))  # PyMOL API
                    cmd.distance("contact_map", "sele1", "sele2")  # PyMOL API
        try:
            cmd.hide("labels", "contact_map")  # PyMOL API
            cmd.delete("sele1")  # PyMOL API
            cmd.delete("sele2")  # PyMOL API
        except:
            pass

    # Show contact map/cross corelation as a graph
    def graph_contact_map(self, plot_type):
        if plot_type == "contact":
            # matplotlib
            prody.showContactMap(self.enm)
        elif plot_type == "cross":
            # matplotlib
            prody.showCrossCorr(self.enm)

    # Show vectors from NMD file
    def show_vectors(self):
        color1 = cmd.get_color_tuple(self.color)  # PyMOL API
        color2 = cmd.get_color_tuple(self.color)  # PyMOL API
        if color1:
            color1 = list(color1)
        # Fallback to grey in case of unrecognized color
        else:
            color1 = [0.5, 0.5, 0.5]
        if color2:
            color2 = list(color2)
        # Fallback to grey in case of unrecognized color
        else:
            color2 = [0.5, 0.5, 0.5]
        arrow_head_radius = 0.15

        x1 = []
        y1 = []
        z1 = []

        coor = "x"
        for coordinate in self.nmd_coordinates:
            if coor == "x":
                x1.append(float(coordinate))
                coor = "y"
            elif coor == "y":
                y1.append(float(coordinate))
                coor = "z"
            elif coor == "z":
                z1.append(float(coordinate))
                coor = "x"

        x2 = []
        y2 = []
        z2 = []

        # This factor is provided to make vector length more like in NMWiz.
        # More investigation is needed to get exact formula.
        approximation_factor = 16.6

        coor = "x"
        coor_nr = 0
        round_nr = 0
        for mode in self.nmd_mode[self.mode_nr]:
            if coor == "x":
                x2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + x1[
                        coor_nr])
                coor = "y"
            elif coor == "y":
                y2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + y1[
                        coor_nr])
                coor = "z"
            elif coor == "z":
                z2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + z1[
                        coor_nr])
                coor = "x"
            round_nr = round_nr + 1
            if round_nr == 3:
                round_nr = 0
                coor_nr = coor_nr + 1

        coor_nr = 0
        for position in x1:
            try:
                cmd.delete("Mode_Vector_" + str(coor_nr))
            except:
                pass

            cone = [cgo.CONE, x1[coor_nr], y1[coor_nr], z1[coor_nr], x2[coor_nr], y2[coor_nr], z2[coor_nr],
                    arrow_head_radius, 0.0] + color1 + color2 + [1.0, 0.0]
            cmd.load_cgo(cone, "Mode_Vector_" + str(coor_nr))  # PyMOL API
            coor_nr = coor_nr + 1
        # Another workaround for PyMOL 1.8 with TravisCI
        try:
            cam_possition = cmd.get_view(quiet=1)  # PyMOL API
            cmd.set_view(cam_possition)  # PyMOL API
        except TypeError:
            pass

    def change_vectors_color(self, color):
        self.color = color
        self.show_vectors()

    def change_vectors_scale(self, scale):
        scale = float(scale)
        self.scale = scale
        self.show_vectors()

    def change_vectors_mode_nr(self, mode_nr):
        self.mode_nr = mode_nr
        self.show_vectors()

    def options_change(self, v1, v2, root):
        self.calculation_type = v1.get()
        self.contact_map = v2.get()
        # save_options()
        root.destroy()

    def block_contact(self, block, contact_map_b, contact_map_v):
        # ToDO: Replace below Tk code with something agnostic
        self.block_contact_map = block
        if block == 0:
            contact_map_b.configure(state=ACTIVE)
        elif block == 1:
            contact_map_b.configure(state=DISABLED)
            contact_map_v.set(0)


# This class create and maintain abstraction mdp file representatives. em.mdp, pr.mdp, md.mdp
class MdpConfig:
    external_file = 0
    options = [[]]
    file_name = ""

    def __init__(self, file_name, init_config, external_file=0):
        self.file_name = file_name
        self.external_file = external_file
        list1 = init_config.split("\n")
        list2 = []
        for line in list1:
            list2.append(line.split(" = "))
        self.options = list2

    def update(self, option_nr, value, check=1):
        self.options[option_nr][1] = value
        if check == 0 and self.options[option_nr][0][0] != ";":
            self.options[option_nr][0] = ";" + self.options[option_nr][0]
        elif check == 1 and self.options[option_nr][0][0] == ";":
            self.options[option_nr][0] = self.options[option_nr][0][1:]
        self.clean_artefacts()

    def save_file(self, s_params):
        project_name = s_params.project_name
        project_dir = get_project_dirs(project_name)
        config = ""
        for option in self.options:
            # pass empty option
            if option == ['']:
                pass
            else:
                config = "{}{} = {}\n".format(config, str(option[0]), str(option[1]))
        mdp = open(project_dir + self.file_name, "w")
        mdp.write(config)
        mdp.close()

    # Clean options from artefacts
    def clean_artefacts(self):
        try:
            self.options.remove([''])
        except:
            pass


# This function creates files needed by the project
def create_config_files(project_name):
    dynamics_dir = get_dynamics_dir()
    project_dir = get_project_dirs(project_name)
    print("Create config files")
    project_dir = get_project_dirs(project_name)
    #    if not os.path.isfile(project_dir + "options.pickle"):
    #        pass
    #    else:
    #        load_options(s_params)
    if os.path.isfile(dynamics_dir + "em.mdp"):
        shutil.copy(dynamics_dir + "em.mdp", project_dir + "em.mdp")
        print("Found em.mdp file. Using it instead of local configuration.")
    elif os.path.isfile(project_dir + "em.mdp"):
        em_file_config = open(project_dir + "em.mdp", "r").read()
        em_file = MdpConfig("em.mdp", em_file_config, 1)
    else:
        em_file = MdpConfig("em.mdp", EM_INIT_CONFIG, 0)
    if os.path.isfile(dynamics_dir + "pr.mdp"):
        shutil.copy(dynamics_dir + "pr.mdp", project_dir + "pr.mdp")
        print("Found pr.mdp file. Using it instead of local configuration.")
    elif os.path.isfile(project_dir + "pr.mdp"):
        pr_file_config = open(project_dir + "pr.mdp", "r").read()
        pr_file = MdpConfig("pr.mdp", pr_file_config, 1)
    else:
        pr_file = MdpConfig("pr.mdp", PR_INIT_CONFIG, 0)
    if os.path.isfile(dynamics_dir + "md.mdp"):
        shutil.copy(dynamics_dir + "md.mdp", project_dir + "md.mdp")
        print("Found md.mdp file. Using it instead of local configuration.")
    elif os.path.isfile(project_dir + "md.mdp"):
        md_file_config = open(project_dir + "md.mdp", "r").read()
        md_file = MdpConfig("md.mdp", md_file_config, 1)
    else:
        md_file = MdpConfig("md.mdp", MD_INIT_CONFIG, 0)
    #    save_options()
    try:
        if project_name in cmd.get_names("objects"):  # PyMOL API
            cmd.save(project_dir + project_name + ".pdb", project_name)  # PyMOL API
            print("cmd saved")
    except (AttributeError, TypeError) as e:
        pass
    return em_file, pr_file, md_file


# Status and to_do maintaining class
class ProgressStatus:
    # 0:Save configuration files; 1:Generate topology file from pdb; 2:Adding Water Box;
    # 3: Adding ions and neutralization 4:Energy Minimization; 5:Position Restrained MD;
    # 6:Restraints; 7:Molecular Dynamics Simulation; 8:Generate multimodel PDB; 9:Calculate vectors using ProDy
    status = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    to_do = [1, 1, 1, 1, 1, 1, 0, 1, 1, 1]

    resume = 0
    x2top = 0
    steps = 8

    def to_do_update(self, position, value):
        if isinstance(position, int):
            self.to_do[position] = value
            self.to_do = self.to_do

    #            save_options()

    def x2top_update(self, value):
        if isinstance(value, int):
            self.x2top = value

    def to_do_status(self):
        to_do = []
        for work in self.status:
            if work == 0:
                to_do.append(1)
            elif work == 1:
                to_do.append(0)
        self.to_do = to_do


# Detect gmx executable along with other associated information...
def get_gromacs_exe_info():
    gmx_exes = ['gmx_mpi_d', 'gmx_mpi', 'gmx']

    gmx_exe = ""
    version = ""
    build_arch = ""
    build_on_cygwin = 0

    stdout_file = "test_gromacs.txt"
    if os.path.isfile(stdout_file):
        os.remove(stdout_file)

    for gmx in gmx_exes:
        cmd = gmx + " -version"
        execute_subprocess(cmd, None, stdout_file)

        ofs = open(stdout_file, "r")
        output = ofs.read()
        ofs.close()

        output = standardize_new_line_char(output)
        if not re.search("GROMACS version:", output, re.I):
            continue

        gmx_exe = gmx
        for line in output.split("\n"):
            if re.search("^[ ]*GROMACS version:", line, re.I):
                gmx_exe = gmx
                version = re.sub("^[ ]*GROMACS version:[ ]*", "", line, flags=re.I)
                if "VERSION " in version:
                    version = version.split("VERSION ")[1].rstrip()
            elif re.search(r"^[ ]*Build OS/arch:", line, re.I):
                build_arch = re.sub("^[ ]*Build OS/arch:[ ]*", "", line, flags=re.I)

                if re.search(r"CYGWIN", build_arch, re.I):
                    build_on_cygwin = 1
        break

    return gmx_exe, version, build_arch, build_on_cygwin


def get_dynamics_dir():
    home_dir = os.path.expanduser('~')
    gmx_home_dir_path = os.path.abspath(home_dir)
    dynamics_dir = os.path.join(gmx_home_dir_path, '.dynamics', '')
    return dynamics_dir


def get_project_dirs(project_name="nothing"):
    dynamics_dir = get_dynamics_dir()
    project_dir = os.path.join(dynamics_dir, project_name, '')
    return project_dir


# Execute command using stdin/stdout as needed...
def execute_subprocess(command, stdin_file_path=None, stdout_file_path=None):
    stdin_file = None
    stdin_msg = "None"
    if stdin_file_path:
        stdin_file = open(stdin_file_path, "r")
        stdin_msg = stdin_file_path

    stdout_file = None
    stdout_msg = "None"
    if stdout_file_path:
        stdout_file = open(stdout_file_path, "w")
        stdout_msg = stdout_file_path

    print("Running command: " + command + "; STDIN: " + stdin_msg + "; STDOUT: " + stdout_msg)

    return_code = subprocess.call(command, stdin=stdin_file, stdout=stdout_file, stderr=subprocess.STDOUT, shell=True)

    if stdin_file_path:
        stdin_file.close()
    if stdout_file_path:
        stdout_file.close()

    return return_code


# Start a subprocess and wait for it to complete along with an option to kill it...
def execute_and_monitor_subprocess(command, stdin_file_path=None, stdout_file_path=None, log_file_path=None):
    if log_file_path:
        if os.path.isfile(log_file_path):
            log_file = open(log_file_path, 'a')
        else:
            log_file = open(log_file_path, 'w')
        star_mark = "\n!{0}!\n".format("*" * 25)
        log_file.write("{0}{1}{0}".format(star_mark, command))
        log_file.close()

    stdin_file = None
    stdin_msg = "None"
    if stdin_file_path:
        stdin_file = open(stdin_file_path, "r")
        stdin_msg = stdin_file_path

    stdout_file = None
    stdout_msg = "None"
    if stdout_file_path:
        stdout_file = open(stdout_file_path, "w")
        stdout_msg = stdout_file_path

    print("Running command: {}; STDIN: {}; STDOUT: {}".format(command, stdin_msg, stdout_msg))

    gmx = subprocess.Popen(command, stdin=stdin_file, stdout=stdout_file, stderr=subprocess.STDOUT, shell=True)
    while gmx.poll() is None:
        #        if stop == 1:
        #            gmx.kill()
        #            break
        time.sleep(1.0)

    if stdin_file_path:
        stdin_file.close()

    if stdout_file_path:
        stdout_file.close()

    # Append any stdout to log file...
    if log_file_path and stdout_file_path:
        log_file = open(log_file_path, "a")
        stdout_file = open(stdout_file_path, "r")
        log_file.write(stdout_file.read())
        log_file.close()
        stdout_file.close()


# Change Windows and Mac new line char to UNIX...
def standardize_new_line_char(in_text):
    out_text = re.sub("(\r\n)|(\r)", "\n", in_text)
    return out_text


# Read text lines and standardize new line char...
def read_text_lines(text_file_path):
    text_lines = []

    ifs = open(text_file_path, "r")
    for line in iter(ifs.readline, ''):
        new_line = standardize_new_line_char(line)
        text_lines.append(new_line)
    ifs.close()

    return text_lines


# Collect water modes information...
def get_water_models_info(gmx_output_lines):
    start_line = 0
    while gmx_output_lines[start_line][0:7] != "Opening":
        start_line = start_line + 1

    start_line = start_line + 1
    end_line = start_line

    while (gmx_output_lines[end_line][0:7] != "Opening") and (gmx_output_lines[end_line][0] != "\n"):
        end_line = end_line + 1

    waters_info = gmx_output_lines[start_line:end_line]
    waters_info2 = []
    number = 1
    for water in waters_info:
        waters_info2.append([number, water[:-1]])
        number = number + 1

    return waters_info2


# Steps mark work as done.
def steps_status_done(step_nr, s_params):
    progress = s_params.progress
    if progress.status[step_nr] == 1:
        return " [done]"
    elif progress.status[step_nr] == 0:
        return ""


# This function will receive status from gromacs2 class and change it to global variable.
def status_update(input_status):
    status = input_status
    print(status[1])


# This function will start real workflow of the plugin, once everything is set
def dynamics(s_params):
    print("Starting PyMOL plugin 'dynamics' ver. {}".format(plugin_ver))
    status = ["ok", ""]
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    progress = s_params.progress
    gromacs2 = s_params.gmx_input
    vectors_prody = s_params.vectors_prody
    os.chdir(project_dir)
    stop = 0

    # Saving configuration files
    if status[0] == "ok" and stop == 0 and progress.to_do[0] == 1:
        mdp_files(s_params)
        if status[0] == "ok":
            progress.status[0] = 1
            progress.to_do[0] = 0
            save_options(s_params)

    # Counting topology
    if status[0] == "ok" and stop == 0 and progress.to_do[1] == 1 and progress.x2top == 0:
        status = gromacs2.pdb2top(s_params)
        if status[0] == "ok":
            progress.status[1] = 1
            progress.to_do[1] = 0
            save_options(s_params)

    elif status[0] == "ok" and stop == 0 and progress.to_do[1] == 1 and progress.x2top == 1:
        status = gromacs2.x2top(s_params)
        if status[0] == "ok":
            progress.status[1] = 1
            progress.to_do[1] = 0
            save_options(s_params)

    # Adding water box
    if status[0] == "ok" and stop == 0 and progress.to_do[2] == 1:
        status = gromacs2.waterbox(s_params)
        if status[0] == "ok":
            progress.status[2] = 1
            progress.to_do[2] = 0
            save_options(s_params)

    # Adding ions
    if status[0] == "ok" and stop == 0 and progress.to_do[3] == 1:
        status = gromacs2.saltadd(s_params)
        if status[0] == "ok":
            progress.status[3] = 1
            progress.to_do[3] = 0
            save_options(s_params)

    # EM
    if status[0] == "ok" and stop == 0 and progress.to_do[4] == 1:
        status = gromacs2.em(s_params)
        if status[0] == "ok":
            progress.status[4] = 1
            progress.to_do[4] = 0
            save_options(s_params)
    elif status[0] == "ok" and stop == 0 and progress.to_do[4] == 0 and progress.status[4] == 0:
        shutil.copy(project_name + "_b4em.gro", project_name + "_b4pr.gro")

    # PR
    if status[0] == "ok" and stop == 0 and progress.to_do[5] == 1:
        status = gromacs2.pr(s_params)
        if status[0] == "ok":
            progress.status[5] = 1
            progress.to_do[5] = 0
            save_options(s_params)
    elif status[0] == "ok" and stop == 0 and progress.to_do[5] == 0 and progress.status[5] == 0:
        shutil.copy(project_name + "_b4pr.gro", project_name + "_b4md.gro")

    # Restraints
    if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
        status = gromacs2.restraints(project_name)
        if status[0] == "ok":
            progress.status[6] = 1
            progress.to_do[6] = 0
            save_options(s_params)

    # MD
    if status[0] == "ok" and stop == 0 and progress.to_do[7] == 1:
        status = gromacs2.md(s_params)
        if status[0] == "ok":
            progress.status[7] = 1
            progress.to_do[7] = 0
            save_options(s_params)

    # Trjconv
    if status[0] == "ok" and stop == 0 and progress.to_do[8] == 1:
        status = gromacs2.trjconv(s_params)
        show_multipdb(s_params)
        progress.status[8] = 1
        progress.to_do[8] = 0
        save_options(s_params)

    # Calculating vectors
    if status[0] == "ok" and stop == 0 and progress.to_do[9] == 1 and prody:
        vectors_prody.prody(project_name)
        vectors_prody.nmd_format(project_name)
        vectors_prody.show_vectors()
        progress.status[9] = 1
        progress.to_do[9] = 0
        save_options(s_params)
    elif status[0] == "fail":
        print(status[1])
        if stop == 0:
            error_message(s_params)


# Saving configuration files
def mdp_files(s_params):
    dynamics_dir = get_dynamics_dir()
    em_file = s_params.em_file
    pr_file = s_params.pr_file
    md_file = s_params.md_file
    if not os.path.isfile("{}em.mdp".format(dynamics_dir)):
        em_file.save_file(s_params)
    if not os.path.isfile("{}pr.mdp".format(dynamics_dir)):
        pr_file.save_file(s_params)
    if not os.path.isfile("{}md.mdp".format(dynamics_dir)):
        md_file.save_file(s_params)


# Show multimodel PDB file in PyMOL
def show_multipdb(s_params):
    project_name = s_params.project_name
    try:
        cmd.hide("everything", project_name)  # PyMOL API
    except (parsing.QuietException, CmdException) as e:  # PyMOL API
        print("Warning: {}".format(e))
    try:
        cmd.load("{}_multimodel.pdb".format(project_name))  # PyMOL API
    except AttributeError:
        pass


# Detect list of PyMOL loaded PDB files if no files than list "nothing"
def get_pdb_names():
    all_names = cmd.get_names("objects")  # PyMOL API
    all_names1 = []
    for name in all_names:
        name1 = name.split("_")
        if name1[-1] == "multimodel" or name1[-1] == "(sele)" or (name1[0] == "Mode" and len(name1) == 3):
            pass
        else:
            all_names1.append(name)
    all_names = all_names1
    if not all_names:
        all_names = ["nothing"]
    return all_names


def read_and_set_init_project(s_params):
    all_names = get_pdb_names()
    project_name = all_names[0]
    s_params.change_project_name(project_name)
    if all_names != ["nothing"]:
        s_params.create_cfg_files()
    return all_names, project_name


# List of previous projects
def list_prev_projects(all_names):
    dynamics_dir = get_dynamics_dir()
    if os.path.isdir(dynamics_dir):
        projects = os.listdir(dynamics_dir)
    else:
        projects = []
    projects2 = []
    for file_dir in projects:
        if os.path.isdir(dynamics_dir + file_dir) and file_dir not in all_names and file_dir != "nothing":
            projects2.append(file_dir)
    return projects2


# Saving tar.bz file
def save_file(destination_path, s_params):
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    print("Saving")
    import tarfile
    save_options(s_params)
    tar = tarfile.open(destination_path + ".tar.bz2", "w:bz2")
    tar.add(project_dir, recursive=True, arcname=project_name)
    tar.close()
    os.remove(destination_path)


# Load tar.bz file
def load_file(file_path, s_params):
    print("Loading file: " + file_path)
    dynamics_dir = get_dynamics_dir()
    tar = tarfile.open(file_path, "r:bz2")
    names = tar.getnames()
    # Backup same name folder if file is loaded
    if os.path.isdir(dynamics_dir + names[0]):
        back_folder = dynamics_dir + names[0] + "_back"
        while os.path.isdir(back_folder):
            back_folder = back_folder + "_b"
        os.rename(dynamics_dir + names[0], back_folder)
    tar.extractall(dynamics_dir)
    project_name = names[0]
    s_params.change_project_name(project_name)
    load_options(s_params)


# Save all settings to options.pickle file
def save_options(s_params):
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    gmx_version = s_params.gmx_output.version
    gromacs2 = s_params.gmx_input
    progress = s_params.progress
    em_file = s_params.em_file
    pr_file = s_params.pr_file
    md_file = s_params.md_file
    if not prody:
        vectors_prody = 0
    else:
        vectors_prody = s_params.vectors_prody
    print("updating project files")
    if not os.path.isdir(project_dir):
        os.makedirs(project_dir)
    destination_option = open(project_dir + "options.pickle", "wb")
    pickle_list = [plugin_ver, gmx_version, gromacs2, em_file, pr_file, md_file, progress, vectors_prody]
    pickle.dump(pickle_list, destination_option)
    del destination_option


# Load all settings from options.pickle file
def load_options(s_params):
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    gmx_version = s_params.gmx_output.version
    gromacs2 = s_params.gmx_input
    pickle_file = open(project_dir + "options.pickle", "rb")
    options = pickle.load(pickle_file)

    print("Loading project {}".format(project_name))
    print("Project was created for Dynamics PyMOL Plugin {} and GROMACS {}".format(options[0], options[1]))
    if gmx_version != options[1]:
        print("GROMACS versions is different for loaded file.")

    if options[0][1:4] == "2.2":
        gromacs2.update({"force": options[2].force, "water": options[2].water, "group": options[2].group,
                         "box_type": options[2].box_type, "hydro": options[2].hydro,
                         "box_distance": options[2].box_distance, "box_density": options[2].box_density,
                         "restraints_nr": options[2].restraints_nr, "neutrality": options[2].neutrality,
                         "salt_conc": options[2].salt_conc, "positive_ion": options[2].positive_ion,
                         "negative_ion": options[2].negative_ion, "explicit": options[2].explicit})
        em_file = options[3]
        pr_file = options[4]
        md_file = options[5]
        progress = options[6]
        if prody and options[7] != 0:
            vectors_prody = options[7]
    elif options[0][1:4] == "2.1":
        print("plugin 2.1 compatibility layer")
        gromacs2.update({"force": options[2].force, "water": options[2].water, "group": options[2].group,
                         "box_type": options[2].box_type, "hydro": options[2].hydro,
                         "box_distance": options[2].box_distance, "box_density": options[2].box_density,
                         "restraints_nr": options[2].restraints_nr, "neutrality": options[2].neutrality,
                         "salt_conc": options[2].salt_conc, "positive_ion": options[2].positive_ion,
                         "negative_ion": options[2].negative_ion})
        em_file = options[3]
        pr_file = options[4]
        md_file = options[5]
        progress = options[6]
        gromacs2.update({"explicit": options[7]})
        if prody and options[8] != 0:
            vectors_prody = options[8]
    else:
        print("Warning. Importing projects from plugin version " + options[0] + " is not supported. Aboring import.")


# Text for "Help"
def help_option():
    help_message = """This is the dynamics PyMOL Plugin.
This software (including its Debian packaging) is available to you under the terms of the GPL-3,
see "/usr/share/common-licenses/GPL-3".
Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
Contributors:
- Tomasz Makarewicz (makson96@gmail.com)
- Ajit B. Datta (ajit@jcbose.ac.in)
- Sara Boch Kminikowska
- Manish Sud (msud@san.rr.com; URL: www.MayaChemTools.org)

Full manual is available to you on project website: https://github.com/makson96/Dynamics/raw/master/manual.odt
or as a file: /usr/share/doc/dynamics-pymol-plugin/manual.odt

The purpose of this plugin is to perform molecular dynamics simulation by GROMACS using easy graphical tool and powerful molecular viewer.

To use this program run it as a PyMOL plugin.
Choose molecule (PDB) for which you want to perform molecular dynamics simulation (left column).
Choose force field and water model options in the middle column.
Choose any additional options in the right column.
Press OK button.
Click Start button and wait till calculation is finished.
Multimodel PDB file will be displayed in PyMOL viewer.
You can click Play button in order to see animation."""
    return help_message


# Clean function
def clean_option():
    shutil.rmtree(get_dynamics_dir())
    print("Temporary files are now removed.")


# If molecular dynamics simulation fails, this function will show the error
def error_message(s_params):
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    log = open(project_dir + "log.txt", "r")
    log_list = log.readlines()
    error_start_line = 0
    for log_line in log_list:
        error_start_line += 1
        if "Fatal error:" in log_line:
            error_start_line -= 0
            break
    error_end_line = error_start_line
    for log_line in log_list[error_start_line:]:
        error_end_line += 1
        if "-------------------------------------------------------" in log_line:
            error_end_line -= 0
            break
    error_list = log_list[error_start_line:error_end_line]
    error = ""
    for line in error_list:
        error = error + line
    print(error)
    # Debug
    file = open("log.txt", "r")
    print(file.read())
    file = open("log1.txt", "r")
    print(file.read())


def create_gui(gui_library, status, s_parameters, parent):
    if status[0] == "ok":
        root_window(status, s_parameters, parent)
    else:
        tkMessageBox.showerror("Initialization error", status[1])

# --Graphic Interface Tk--
# Don't care too much of below code quality, as Tk it depreciated and will be removed in plugin version 3.1
# Root menu window
def root_window(status, s_params, parent):
    if parent:
        root = parent
    else:
        # First try to get this root fails, but the second try works fine.
        try:
            root_pymol = plugins.get_tk_root()
        except ModuleNotFoundError:
            root_pymol = plugins.get_tk_root()
        root = Toplevel(root_pymol)
    root.wm_title("Dynamics with Gromacs" + plugin_ver)
    calculationW = CalculationWindow()
    waterW = WaterWindows()
    restraintsW = RestraintsWindow()
    genionW = GenionWindow()
    gromacs = s_params.gmx_output
    gromacs2 = s_params.gmx_input
    dynamics_dir = get_dynamics_dir()
    vectors_prody = s_params.vectors_prody
    all_names, project_name = read_and_set_init_project(s_params)

    # TkInter variables
    v1_name = StringVar(root)
    v1_name.set(project_name)

    group_nr = gromacs.group_list[1][0]
    v2_group = IntVar(root)
    v2_group.set(group_nr)

    force_nr = gromacs.force_list[0][0]
    v3_force = IntVar(root)
    v3_force.set(force_nr)

    water_nr = gromacs.water_list[0][0]
    v4_water = IntVar(root)
    v4_water.set(water_nr)

    water_v = StringVar(root)
    water_v.set(gromacs.water_list[0][1])

    time_entry_value = StringVar(root)
    time_entry_value.set("10.0")

    # Start drawing interface
    frame0 = Frame(root)
    frame0.pack(side=TOP)

    w_version = Label(frame0, text="GROMACS VERSION " + gromacs.version)
    w_version.pack(side=TOP)

    frame1 = Frame(root)
    frame1.pack(side=TOP)

    frame1_1 = Frame(frame1, borderwidth=1, relief=RAISED)
    frame1_1.pack(side=LEFT)

    w1 = Label(frame1_1, text="Molecules", font="bold")
    w1.pack(side=TOP)

    frame1_1a = Frame(frame1_1)
    frame1_1a.pack(side=TOP)

    # List of PyMOL loaded PDB files
    if all_names[0] != "nothing":
        for molecule in all_names:
            radio_button1 = Radiobutton(frame1_1a, text=molecule, value=molecule, variable=v1_name,
                                        command=lambda: set_variables(v1_name.get(), v2_group, v3_force, v4_water,
                                                                      water_v, check1_button, s_params))
            radio_button1.pack(side=TOP, anchor=W)
    # If no loaded PDB files, than add button to choose one
    else:
        w1_1 = Label(frame1_1a, text="Choose PDB file")
        w1_1.pack(side=TOP)
        frame1_1_1 = Frame(frame1_1a)
        frame1_1_1.pack(side=TOP)
        label1 = Label(frame1_1_1, textvariable=v1_name)
        label1.pack(side=LEFT)
        button_e1 = Button(frame1_1_1, text="Browse", command=lambda: select_file(v1_name, s_params))
        button_e1.pack(side=LEFT)

    # List of previous projects
    projects = list_prev_projects(all_names)
    if projects:
        w1_2 = Label(frame1_1, text="Previous Projects")
        w1_2.pack(side=TOP)

        for molecule in projects:
            molecule1 = molecule.split("_")
            if molecule1[-1] == "multimodel":
                pass
            else:
                molecule1 = molecule.split("-")
                if molecule1[0] == "gromacs":
                    pass
                else:
                    radio_button1 = Radiobutton(frame1_1, text=molecule, value=molecule, variable=v1_name,
                                                command=lambda: set_variables(v1_name.get(), v2_group, v3_force,
                                                                              v4_water, water_v, check1_button))
                    radio_button1.pack(side=TOP, anchor=W)

    # List of group for final model
    w2 = Label(frame1_1, text="Group", font="bold")
    w2.pack(side=TOP)
    for group in gromacs.group_list:
        radio_button2 = Radiobutton(frame1_1, text=group[1], value=group[0], variable=v2_group,
                                    command=lambda: gromacs2.update({"group": v2_group.get()}))
        radio_button2.pack(side=TOP, anchor=W)

    frame1_2 = Frame(frame1, borderwidth=1, relief=RAISED)
    frame1_2.pack(side=LEFT)

    # List of available force fields
    w3 = Label(frame1_2, text="Force fields", anchor=E, font="bold")
    w3.pack(side=TOP)

    for force in gromacs.force_list:
        radio_button3 = Radiobutton(frame1_2, text=force[1], value=force[0], variable=v3_force,
                                    command=lambda: waterW.change(v4_water, water_v, v3_force.get()))
        radio_button3.pack(side=TOP, anchor=W)

    # Label of choosen water model
    w4 = Label(frame1_2, text="Water Model", anchor=E, font="bold")
    w4.pack(side=TOP)

    frame1_2_1 = Frame(frame1_2)
    frame1_2_1.pack(side=TOP)

    # Buttons to choose water model and configure water box
    water_label = Label(frame1_2_1, textvariable=water_v)
    water_label.pack(side=LEFT)
    water_button = Button(frame1_2_1, text="Choose...",
                          command=lambda: waterW.choose(v4_water, water_v, waterbox_button, root, s_params))
    water_button.pack(side=LEFT)
    waterbox_button = Button(frame1_2_1, text="Configure", command=lambda: waterW.box(root, s_params))
    waterbox_button.pack(side=LEFT)
    waterbox_button2 = Button(frame1_2_1, text="Hydrogen Mass", command=lambda: waterW.box2(root, s_params))
    waterbox_button2.pack(side=LEFT)

    frame1_3 = Frame(frame1)
    frame1_3.pack(side=LEFT)

    frame1_3_1 = Frame(frame1_3, borderwidth=1, relief=RAISED)
    frame1_3_1.pack(side=TOP)

    w4 = Label(frame1_3_1, text="Configuration", font="bold")
    w4.pack(side=TOP)

    # Button for configuration of Simulation Steps
    steps_label = Label(frame1_3_1, text="Simulation Steps")
    steps_label.pack(side=TOP)
    steps_button = Button(frame1_3_1, text="Configure",
                          command=lambda: steps_configure(root, check1_button, s_params, restraintsW))
    steps_button.pack(side=TOP)

    # Button for Genion configuration
    ion_label = Label(frame1_3_1, text="Adding ions & Neutralize")
    ion_label.pack(side=TOP)
    ion_button2 = Button(frame1_3_1, text="Configure", command=lambda: genionW.window(root, s_params))
    ion_button2.pack(side=TOP)

    # Button for configuration of MDP files
    em_label = Label(frame1_3_1, text="Energy Minimization")
    em_label.pack(side=TOP)
    em_button2 = Button(frame1_3_1, text="Configure", command=lambda: mdp_configure("em", root, s_params))
    em_button2.pack(side=TOP)
    if os.path.isfile(dynamics_dir + "em.mdp"):
        em_button2.configure(state=DISABLED)

    pr_label = Label(frame1_3_1, text="Position Restrained MD")
    pr_label.pack(side=TOP)
    pr_button2 = Button(frame1_3_1, text="Configure", command=lambda: mdp_configure("pr", root, s_params))
    pr_button2.pack(side=TOP)
    if os.path.isfile(dynamics_dir + "pr.mdp"):
        pr_button2.configure(state=DISABLED)

    md_label = Label(frame1_3_1, text="Molecular Dynamics Simulation")
    md_label.pack(side=TOP)
    md_button2 = Button(frame1_3_1, text="Configure", command=lambda: mdp_configure("md", root, s_params))
    md_button2.pack(side=TOP)
    if os.path.isfile(dynamics_dir + "md.mdp"):
        md_button2.configure(state=DISABLED)

    # Button for configuration of Restraints
    re_label = Label(frame1_3_1, text="Restraints (Select Atoms)")
    re_label.pack(side=TOP)

    check1_button = Button(frame1_3_1, text="Configure", command=lambda: restraintsW.window(root, s_params))
    check1_button.pack(side=TOP)

    # Button for ProDy options
    pro_label = Label(frame1_3_1, text="Vectors Options")
    pro_label.pack(side=TOP)
    prody_button = Button(frame1_3_1, text="Configure", command=lambda: vectors_prody.window(root))
    prody_button.pack(side=TOP)

    # Dynamics Simulation Time
    time_label = Label(frame1_3_1, text="Dynamics Simulation Time")
    time_label.pack(side=TOP)

    frame1_3_1_1 = Frame(frame1_3_1)
    frame1_3_1_1.pack(side=TOP)

    time_entry = Entry(frame1_3_1_1, textvariable=time_entry_value)
    time_entry.pack(side=LEFT)
    time_label2 = Label(frame1_3_1_1, text="[ps]")
    time_label2.pack(side=LEFT)
    time_button = Button(frame1_3_1_1, text="OK", command=lambda: s_params.md_file.update(3, str(
        int(float(time_entry_value.get()) / float(s_params.md_file.options[2][1])))))
    time_button.pack(side=LEFT)

    # Disable configuration of ProDy (Vectors) if ProDy is not installed
    if not prody:
        prody_button.configure(state=DISABLED)

    frame2 = Frame(root)
    frame2.pack(side=TOP)

    # Additional Buttons
    exit_button = Button(frame2, text="Exit", command=root.destroy)
    exit_button.pack(side=LEFT)

    clean_button = Button(frame2, text="Clean", command=clean_message)
    clean_button.pack(side=LEFT)

    help_button = Button(frame2, text="Help", command=lambda: help_window(root))
    help_button.pack(side=LEFT)

    save_button = Button(frame2, text="Save", command=select_file_save)
    save_button.pack(side=LEFT)

    load_button = Button(frame2, text="Load",
                         command=lambda: select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v,
                                                          check1_button, s_params))
    load_button.pack(side=LEFT)

    count_button = Button(frame2, text="OK", command=lambda: calculationW.check_window(root, root_pymol, s_params,
                                                                                       status))
    count_button.pack(side=LEFT)

    # Initial configuration
    set_variables(v1_name.get(), v2_group, v3_force, v4_water, water_v, check1_button, s_params)


# Molecular Dynamics Performing window
class CalculationWindow:
    tasks_to_do = 0
    bar_var = ""
    bar_widget = ""
    start_button = ""
    stop_button = ""
    log_button = ""

    def __init__(self):
        self.queue_status = Queue.Queue()
        self.queue_percent = Queue.Queue()

    # This will prevent Calculation Window to display if non protein has been selected
    def check_window(self, master, g_parent, s_params, status):
        project_name = s_params.project_name
        if project_name != "nothing":
            master.destroy()
            root = Toplevel(g_parent)
            self.window(root, s_params, status, g_parent)
        elif project_name == "nothing":
            no_molecule_warning()

    # This function will create main Calculation Window
    def window(self, root, s_params, status, parent):
        root.wm_title("Calculation Window")
        frame1 = Frame(root)
        frame1.pack(side=TOP)
        frame2 = Frame(root)
        frame2.pack(side=TOP)

        self.bar_var = StringVar(root)
        self.bar_var.set("Ready to start")

        w5 = Label(frame1, textvariable=self.bar_var)
        w5.pack(side=TOP)
        self.bar_widget = Progressbar(frame1)
        self.bar_widget.pack(side=TOP)

        exit_button = Button(frame2, text="EXIT", command=root.destroy)
        exit_button.pack(side=LEFT)

        save_button = Button(frame2, text="SAVE", command=lambda: select_file_save(1))
        save_button.pack(side=LEFT)

        stop_button = Button(frame2, text="STOP", command=lambda: self.start_counting(0))
        stop_button.pack(side=LEFT)
        stop = s_params.stop
        if stop:
            stop_button.configure(state=DISABLED)
        self.stop_button = stop_button

        start_button = Button(frame2, text="START", command=lambda: self.start_counting(1))
        start_button.pack(side=LEFT)
        if stop == 0:
            start_button.configure(state=DISABLED)
        self.start_button = start_button

        log_button = Button(frame2, text="LOG", command=log_window)
        log_button.pack(side=LEFT)
        log_button.configure(state=DISABLED)
        self.log_button = log_button

        # Updateing status bar
        tasks_nr = 0.0
        for task in s_params.progress.to_do:
            tasks_nr = tasks_nr + task
        self.tasks_to_do = tasks_nr
        thread.start_new_thread(self.bar_update, s_params, status)
        self.bar_display(root, parent, s_params)

    # This function will update status bar during molecular dynamics simulation (beware this is separate thread)
    def bar_update(self, s_params, status):
        percent = 0.0
        while s_params.stop:
            time.sleep(0.5)
        while percent != 100:  # and error == ""
            time.sleep(0.5)
            percent = steps_status_bar("only_bar", s_params)
            self.queue_percent.put(percent)
            if s_params.stop == 0:
                self.queue_status.put(status[1])
            elif s_params.stop == 1:
                self.queue_status.put("User Stoped")

    #        if error != "":
    #            self.queue_status.put("Fatal Error")

    # This function will update status bar in thread safe manner
    def bar_display(self, root, parent, s_params):
        try:
            status = self.queue_status.get(block=False)
            self.bar_var.set(status)
        except Queue.Empty:
            status = "No change"
        try:
            percent = self.queue_percent.get(block=False)
            self.bar_widget.configure(value=percent)
        except:
            pass
        if status == "Fatal Error":
            self.start_counting(0)
            self.start_button.configure(state=DISABLED)
            tkMessageBox.showerror("GROMACS Error Message", "Error")  # error)
        if status == "Finished!":
            root.destroy()
            # Show interpretation window after successful completion of the calculations...
            show_interpretation_window(parent, s_params)
        else:
            root.after(100, self.bar_display, root)

    # This function will change global value if stop is clicked during simulation
    def start_counting(self, value):
        if value == 1:
            stop = 0
            thread.start_new_thread(dynamics, ())
            self.stop_button.configure(state=ACTIVE)
            self.start_button.configure(state=DISABLED)
            self.log_button.configure(state=DISABLED)
        elif value == 0:
            stop = 1
            self.stop_button.configure(state=DISABLED)
            self.start_button.configure(state=ACTIVE)
            self.log_button.configure(state=ACTIVE)


# This window will allow to manipulate final molecule to interprate MD simulation results
class InterpretationWindow:
    dt = 0.0
    nsteps = 0.0
    nstxout = 0.0
    max_time = 0.0
    tentry_value = ""
    pause = 1

    def __init__(self, g_parent, s_params):
        self.queue_time = Queue.Queue()
        self.md_time(s_params)
        root = Toplevel(g_parent)
        self.window(root, s_params)

    def md_time(self, s_params):
        project_name = s_params.project_name
        project_dir = get_project_dirs(project_name)
        md_file = open(project_dir + "md.mdp", "r")
        for lines in md_file.readlines():
            splited_line = lines.split(" ")
            if splited_line[0] == "dt":
                dt = float(splited_line[2])
                self.dt = dt
            elif splited_line[0] == "nsteps":
                nsteps = float(splited_line[2])
                self.nsteps = nsteps
            elif splited_line[0] == "nstxout":
                nstxout = float(splited_line[2])
                self.nstxout = nstxout
        max_time = dt * nsteps
        self.max_time = max_time

    def window(self, root, s_params):
        vectors_prody = s_params.vector_prody
        root.wm_title("MD Interpretation")
        self.tentry_value = StringVar(root)
        self.tentry_value.set("0.0")
        sentry_value = StringVar(root)
        sentry_value.set("1.0")
        contact_entry_value = StringVar(root)
        contact_entry_value.set("-1.0")

        frame1 = Frame(root)
        frame1.pack(side=TOP)

        # Animation
        alabel = Label(frame1, text="Animation", font="bold")
        alabel.pack()
        frame1_1 = Frame(frame1)
        frame1_1.pack(side=TOP)
        play_button = Button(frame1_1, text="PLAY", command=lambda: self.pause_play(0))
        play_button.pack(side=LEFT)
        pause_button = Button(frame1_1, text="PAUSE", command=lambda: self.pause_play(1))
        pause_button.pack(side=LEFT)
        frame1_2 = Frame(frame1)
        frame1_2.pack(side=TOP, anchor=W)
        tlabel = Label(frame1_2, text="Time [ps] (Max " + str(self.max_time) + " [ps])")
        tlabel.pack(side=LEFT)
        tentry = Entry(frame1_2, textvariable=self.tentry_value)
        tentry.pack(side=LEFT)
        tok_button = Button(frame1_2, text="OK",
                            command=lambda: cmd.frame(self.time2frames(self.tentry_value.get())))  # PyMOL API
        tok_button.pack(side=LEFT)
        frame1_3 = Frame(frame1)
        frame1_3.pack(side=TOP, anchor=W)
        mlabel = Label(frame1_3, text="Model Type")
        mlabel.pack(side=LEFT)
        lines_button = Button(frame1_3, text="Lines", command=lambda: self.shape("lines", s_params))
        lines_button.pack(side=LEFT)
        sticks_button = Button(frame1_3, text="Sticks", command=lambda: self.shape("sticks", s_params))
        sticks_button.pack(side=LEFT)
        ribbon_button = Button(frame1_3, text="Ribbon", command=lambda: self.shape("ribbon", s_params))
        ribbon_button.pack(side=LEFT)
        cartoon_button = Button(frame1_3, text="Cartoon", command=lambda: self.shape("cartoon", s_params))
        cartoon_button.pack(side=LEFT)
        frame1_3_1 = Frame(frame1)
        frame1_3_1.pack(side=TOP, anchor=W)
        mlabel = Label(frame1_3_1, text="Labels")
        mlabel.pack(side=LEFT)
        end_button = Button(frame1_3_1, text="Terminus", command=lambda: self.label("terminus", s_params))
        end_button.pack(side=LEFT)
        acids_button = Button(frame1_3_1, text="Amino Acids", command=lambda: self.label("acids", s_params))
        acids_button.pack(side=LEFT)
        clear_button = Button(frame1_3_1, text="Clear", command=lambda: self.label("clear", s_params))
        clear_button.pack(side=LEFT)

        thread.start_new_thread(self.watch_frames, ())
        self.display_time(root)

        # Vectors
        vlabel = Label(frame1, text="Vectors (Require ProDy)", font="bold")
        vlabel.pack()
        frame1_4 = Frame(frame1)
        frame1_4.pack(side=TOP, anchor=W)
        modlabel = Label(frame1_4, text="Mode Nr")
        modlabel.pack(side=LEFT)
        one_button = Button(frame1_4, text="1", command=lambda: vectors_prody.change_vectors_mode_nr(0))
        one_button.pack(side=LEFT)
        two_button = Button(frame1_4, text="2", command=lambda: vectors_prody.change_vectors_mode_nr(1))
        two_button.pack(side=LEFT)
        three_button = Button(frame1_4, text="3", command=lambda: vectors_prody.change_vectors_mode_nr(2))
        three_button.pack(side=LEFT)
        frame1_5 = Frame(frame1)
        frame1_5.pack(side=TOP, anchor=W)
        slabel = Label(frame1_5, text="Scale")
        slabel.pack(side=LEFT)
        sentry = Entry(frame1_5, textvariable=sentry_value)
        sentry.pack(side=LEFT)
        sok_button = Button(frame1_5, text="OK", command=lambda: vectors_prody.change_vectors_scale(sentry_value.get()))
        sok_button.pack(side=LEFT)
        frame1_6 = Frame(frame1)
        frame1_6.pack(side=TOP, anchor=W)
        modlabel = Label(frame1_6, text="Color")
        modlabel.pack(side=LEFT)
        gray_button = Button(frame1_6, text="Gray", command=lambda: vectors_prody.change_vectors_color("gray"))
        gray_button.pack(side=LEFT)
        red_button = Button(frame1_6, text="Red", command=lambda: vectors_prody.change_vectors_color("red"))
        red_button.pack(side=LEFT)
        blue_button = Button(frame1_6, text="Blue", command=lambda: vectors_prody.change_vectors_color("blue"))
        blue_button.pack(side=LEFT)
        green_button = Button(frame1_6, text="Green", command=lambda: vectors_prody.change_vectors_color("green"))
        green_button.pack(side=LEFT)

        frame1_7 = Frame(frame1)
        frame1_7.pack(side=TOP, anchor=W)
        modlabel = Label(frame1_7, text="Plot results")
        modlabel.pack(side=LEFT)
        contact_button = Button(frame1_7, text="Show Contact Map Graph",
                                command=lambda: vectors_prody.graph_contact_map("contact"))
        contact_button.pack(side=LEFT)
        cross_button = Button(frame1_7, text="Show Cross-correlations Graph",
                              command=lambda: vectors_prody.graph_contact_map("cross"))
        cross_button.pack(side=LEFT)
        frame1_8 = Frame(frame1)
        frame1_8.pack(side=TOP, anchor=W)
        modlabel = Label(frame1_8, text="Plot results")
        modlabel.pack(side=LEFT)
        contact_pymol_button = Button(frame1_8, text="Show Contact Map In PyMOL",
                                      command=lambda: vectors_prody.show_contact_map(contact_entry_value.get(),
                                                                                     s_params.project_name))
        contact_pymol_button.pack(side=LEFT)
        contact_label = Label(frame1_8, text="Sensitivity")
        contact_label.pack(side=LEFT)
        contact_entry = Entry(frame1_8, textvariable=contact_entry_value)
        contact_entry.pack(side=LEFT)

        frame1_8 = Frame(frame1)
        frame1_8.pack(side=TOP)
        exit_button = Button(frame1_8, text="Exit", command=root.destroy)
        exit_button.pack(side=LEFT)
        save_button = Button(frame1_8, text="Save", command=lambda: select_file_save(s_params))
        save_button.pack(side=LEFT)
        log_button = Button(frame1_8, text="Log", command=log_window)
        log_button.pack(side=LEFT)

        if not prody:
            print("No ProDy found")
            one_button.configure(state=DISABLED)
            two_button.configure(state=DISABLED)
            three_button.configure(state=DISABLED)
            sok_button.configure(state=DISABLED)
            gray_button.configure(state=DISABLED)
            red_button.configure(state=DISABLED)
            blue_button.configure(state=DISABLED)
            green_button.configure(state=DISABLED)
        if not prody or vectors_prody.contact_map != 1:
            contact_button.configure(state=DISABLED)
            cross_button.configure(state=DISABLED)
            contact_pymol_button.configure(state=DISABLED)

    def pause_play(self, value):
        if value == 1:
            self.pause = 1
            cmd.mstop()  # PyMOL API
        elif value == 0:
            self.pause = 0
            cmd.mplay()  # PyMOL API

    def frames2time(self, text_var):
        frame = float(text_var)
        time = frame * self.dt * self.nstxout
        return time

    def time2frames(self, text_var):
        nsecond = float(text_var)
        frame = nsecond / self.dt / self.nstxout
        frame = int(frame)
        return frame

    @staticmethod
    def shape(shape_type, s_params):
        project_name = s_params.project_name
        cmd.hide("everything", project_name + "_multimodel")  # PyMOL API
        cmd.show(shape_type, project_name + "_multimodel")  # PyMOL API

    @staticmethod
    def label(name, s_params):
        project_name = s_params.project_name
        if name == "terminus":
            cmd.label("n. ca and {}_multimodel and i. 1".format(project_name), '"N-terminus"')  # PyMOL API
            ca_number = cmd.count_atoms("n. ca and " + project_name + "_multimodel")  # PyMOL API
            cmd.label("n. ca and {}_multimodel and i. {}".format(project_name, str(ca_number)),
                      '"C-terminus"')  # PyMOL API
        elif name == "acids":
            cmd.label("n. ca and {}_multimodel".format(project_name), "resn")  # PyMOL API
        elif name == "clear":
            cmd.label("n. ca and {}_multimodel".format(project_name), "")  # PyMOL API

    # This function will watch time (beware this is separate thread)
    def watch_frames(self):
        while 1:
            pymol_frame = cmd.get_frame()  # PyMOL API
            pymol_time = self.frames2time(pymol_frame)
            self.queue_time.put(pymol_time)
            time.sleep(0.1)

    # This function will update display time in thread safe manner
    def display_time(self, root):
        try:
            time = self.queue_time.get(block=False)
        except Queue.Empty:
            time = "No change"
        if self.pause != 1:
            self.tentry_value.set(time)
        root.after(100, self.display_time, root)


# Show interpretation window...
def show_interpretation_window(parent, s_params):
    InterpretationWindow(parent, s_params)


def help_window(master):
    root = Toplevel(master)
    root.wm_title("Help Window")
    frame = Frame(root)
    frame.pack()

    w = Label(frame, text=help_option())
    w.pack()
    ok_button = Button(frame, text="OK", command=root.destroy)
    ok_button.pack()


def log_window(s_params):
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    if sys.platform == "linux2":
        cmd = "xdg-open {}log.txt".format(project_dir)
        execute_subprocess(cmd)
    elif sys.platform == "darwin":
        cmd = "open {}log.txt".format(project_dir)
        execute_subprocess(cmd)
    elif sys.platform.startswith('win'):
        cmd = "start {}log.txt".format(project_dir)
        execute_subprocess(cmd)


def clean_message():
    tkMessageBox.showinfo("Clean", "Temporary files are now removed!\nPlease restart plugin.")
    clean_option()


def no_molecule_warning():
    tkMessageBox.showinfo("No Molecule Selected", "Please choose any molecule before using this option.")


# This class is resposible for graphic edition of restraints
class RestraintsWindow:
    atom_list = []
    check_var = ""

    # This function will create main window for restraints
    def window(self, master, s_params):
        gromacs2 = s_params.gmx_input
        gromacs = s_params.gmx_output
        root = Toplevel(master)
        root.wm_title("Restraints Configure")

        ok_button = Button(root, text="OK", command=lambda: self.index(s_params))
        ok_button.pack(side=BOTTOM)

        sb = Scrollbar(root, orient=VERTICAL)
        sb.pack(side=RIGHT, fill=Y)

        canvas = Canvas(root, width=600)
        canvas.pack(side=TOP, fill="both", expand=True)
        frame1 = Frame(canvas)
        frame1.pack(side=TOP)

        # attach canvas (with frame1 in it) to scrollbar
        canvas.config(yscrollcommand=sb.set)
        sb.config(command=canvas.yview)

        # bind frame1 with canvas
        canvas.create_window((1, 1), window=frame1, anchor="nw", tags="frame1")
        frame1.bind("<Configure>", canvas.config(scrollregion=(0, 0, 0, 4500)))

        self.check_var = IntVar(frame1)
        self.check_var.set(gromacs2.restraints_nr)

        self.atom_list = []
        number = 0
        for group in gromacs.restraints:
            select = Radiobutton(frame1, text=group[0], value=number, variable=self.check_var)
            select.pack()
            text = Text(frame1)
            text.insert(END, group[1])
            text.pack()
            self.atom_list.append(text)
            number = number + 1

        select1 = Radiobutton(frame1, text="[ PyMol Selected ]", value=number, variable=self.check_var)
        select1.pack()
        text1 = Text(frame1)

        stored.list = []
        cmd.iterate("(sele)", "stored.list.append(ID)")  # PyMOL API

        stored_string = ""
        for atom in stored.list:
            stored_string = stored_string + str(atom)
            lengh = stored_string.split('\n')
            if len(lengh[-1]) < 72:
                stored_string = stored_string + "   "
            else:
                stored_string = stored_string + "\n"

        text1.insert(END, stored_string)
        text1.pack()
        self.atom_list.append(text1)

    # This function will modyfie index_dynamics.ndx file based on user choosed restraints
    def index(self, s_params, root_to_kill=False):
        gromacs2 = s_params.gmx_input
        gromacs = s_params.gmx_output
        index_nr = self.check_var.get()
        gromacs2.restraints_nr = index_nr
        text = self.atom_list[index_nr]
        if index_nr < len(gromacs.restraints):
            gromacs.restraints[index_nr][1] = text.get(1.0, END)
            gromacs.restraints = gromacs.restraints
        index_file = open("index_dynamics.ndx", "w")
        index_file.write("[ Dynamics Selected ]\n" + text.get(1.0, END))
        index_file.close()
        if root_to_kill:
            root_to_kill.destroy()

    # This function will activ or disable restraints button in main window based on check box
    def check(self, check, config_button, s_params):
        md_file = s_params.md_file
        gromacs = s_params.gmx_output
        progress = s_params.progress
        if check == 1:
            config_button.configure(state=ACTIVE)
            md_file.update(2, md_file.options[2][1], 1)
            gromacs.restraints_index()
            progress.to_do[6] = 1
            progress.to_do = progress.to_do
        elif check == 0:
            config_button.configure(state=DISABLED)
            md_file.update(2, md_file.options[2][1], 0)
            progress.to_do[6] = 0
            progress.to_do = progress.to_do


# This function will create window, which allow you to choose PDB file if no file is loaded to PyMOL
def select_file(v_name, s_params):
    root = Tk()
    file = tkFileDialog.askopenfile(parent=root, mode='rb', title='Choose PDB file')
    try:
        name = file.name.split("/")
        name2 = name[-1].split(".")
        # Checking directories
        project_name = name2[0]
        s_params.change_project_name(project_name)
        project_dir = get_project_dirs(project_name)
        v_name.set(project_name)
        if not os.path.isdir(project_dir):
            os.makedirs(project_dir)
            shutil.copyfile(file.name, project_dir + project_name + ".pdb")
            print("pdb_copied")
        create_config_files(project_name)
    except:
        pass
    root.destroy()


# This function will create window, which allow you to save current work
def select_file_save(s_params, rest_of_work=0):
    project_name = s_params.project_name
    progress = s_params.progress
    if project_name != "nothing":
        if rest_of_work == 1:
            progress.to_do_status()
        root = Tk()
        file = tkFileDialog.asksaveasfile(parent=root, mode='w', title='Choose save file')
        if not file:
            save_file(file.name, s_params)
        root.destroy()
    elif project_name == "nothing":
        no_molecule_warning()


# This function will create window, which allow you to load previously saved work
def select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v, config_button_restraints, s_params):
    project_name = s_params.project_name
    gromacs = s_params.gmx_output
    gromacs2 = s_params.gmx_input
    root = Tk()
    file = tkFileDialog.askopenfile(parent=root, mode='rb', defaultextension=".tar.bz2", title='Choose file to load')
    if not file:
        load_file(file.name, s_params)
        v1_name.set(project_name)
        v2_group.set(gromacs.group_list[gromacs2.group][0])
        v3_force.set(gromacs.force_list[gromacs2.force - 1][0])
        v4_water.set(gromacs.water_list[gromacs2.water - 1][0])
        water_v.set(gromacs.water_list[v4_water.get() - 1][1])
        radio_button1 = Radiobutton(frame1_1a, text=project_name, value=project_name, variable=v1_name,
                                    command=lambda: set_variables(v1_name.get(), v2_group, v3_force, v4_water, water_v,
                                                                  config_button_restraints))
        radio_button1.pack(side=TOP, anchor=W)
    root.destroy()


# This function sets variables after choosing new molecule
def set_variables(name, v2_group, v3_force, v4_water, water_v, config_button_restraints, s_params):
    print("Set Variables")
    gromacs = s_params.gmx_output
    gromacs2 = s_params.gmx_input
    progress = s_params.progress
    # Set project name and dir
    project_name = name
    if name:
        s_params.change_project_name(project_name)
    project_dir = get_project_dirs(project_name)
    if os.path.isfile("{}options.pickle".format(project_dir)):
        load_options(s_params)
        v2_group.set(gromacs.group_list[gromacs2.group][0])
        v3_force.set(gromacs.force_list[gromacs2.force - 1][0])
        v4_water.set(gromacs.water_list[gromacs2.water - 1][0])
        water_v.set(gromacs.water_list[v4_water.get() - 1][1])
    else:
        create_config_files(project_name)
    # Correct set of restraints button
    if progress.to_do[6] == 0:
        config_button_restraints.configure(state=DISABLED)
    elif progress.to_do[6] == 1:
        config_button_restraints.configure(state=ACTIVE)
    # If Resume is zero than initial Steps are all ON
    if progress.resume == 0:
        progress.to_do = [1, 1, 1, 1, 1, 1, 0, 1, 1, 1]


# This function will create the window with configuration files based on MDP class
def mdp_configure(config_name, master, s_params):
    project_name = s_params.project_name
    em_file = s_params.em_file
    pr_file = s_params.pr_file
    md_file = s_params.md_file
    if project_name != "nothing":
        root2 = Toplevel(master)

        if config_name == "em":
            em_file.clean_artefacts()
            options = em_file.options
            root2.wm_title("Energy Minimization Options")
        elif config_name == "pr":
            pr_file.clean_artefacts()
            options = pr_file.options
            root2.wm_title("Position Restrained MD Options")
        elif config_name == "md":
            md_file.clean_artefacts()
            options = md_file.options
            root2.wm_title("Molecular Dynamics Simulation Options")

        values_list = []
        check_list = []

        if config_name == "em":
            b = Button(root2, text="OK", command=lambda: mdp_update(values_list, check_list, "em", s_params, root2))
            b.pack(side=BOTTOM)
        elif config_name == "pr":
            b = Button(root2, text="OK", command=lambda: mdp_update(values_list, check_list, "pr", s_params, root2))
            b.pack(side=BOTTOM)
        elif config_name == "md":
            b = Button(root2, text="OK", command=lambda: mdp_update(values_list, check_list, "md", s_params, root2))
            b.pack(side=BOTTOM)

        sb = Scrollbar(root2, orient=VERTICAL)
        sb.pack(side=RIGHT, fill=Y)

        canvas = Canvas(root2, width=400)
        canvas.pack(side=TOP, fill="both", expand=True)
        frame1 = Frame(canvas)
        frame1.pack(side=TOP)

        # attach canvas (with frame1 in it) to scrollbar
        canvas.config(yscrollcommand=sb.set)
        sb.config(command=canvas.yview)

        # bind canvas with frame1 1/2
        canvas.create_window((1, 1), window=frame1, anchor="nw", tags="frame1")
        for option, value in options:
            frame2 = Frame(frame1)
            frame2.pack(side=TOP)
            if option == "emtol":
                l1 = Label(frame2, text="Energy minimizing stuff")
                l1.pack(side=TOP)
            elif option == "Tcoupl":
                l1 = Label(frame2, text="Berendsen temperature and coupling")
                l1.pack(side=TOP)
            elif option == "Pcoupl":
                l1 = Label(frame2, text="Pressure coupling")
                l1.pack(side=TOP)
            elif option == "gen_vel":
                l1 = Label(frame2, text="Generate velocites temperature")
                l1.pack(side=TOP)
            elif option == "constraints":
                l1 = Label(frame2, text="Options for bonds")
                l1.pack(side=TOP)
            values_list.append(StringVar(root2))
            values_list[-1].set(value)
            check_list.append(IntVar(root2))
            if option[0] != ";":
                check_list[-1].set(1)
                c1 = Checkbutton(frame2, text=option, variable=check_list[-1], width=25, anchor=W)
                c1.pack(side=LEFT)
            else:
                check_list[-1].set(0)
                c1 = Checkbutton(frame2, text=option, variable=check_list[-1], width=25, anchor=W)
                c1.pack(side=LEFT)
            e = Entry(frame2, textvariable=values_list[-1])
            e.pack(side=LEFT)

        # bind canvas with frame1 2/2
        frame1.bind("<Configure>", canvas.config(scrollregion=(0, 0, 0, len(values_list) * 25)))

    elif project_name == "nothing":
        no_molecule_warning()


# This function will update MDP class objects alfter closing "mdp_configure" window
def mdp_update(values, check_list, mdp, s_params, root_to_kill=""):
    em_file = s_params.em_file
    pr_file = s_params.pr_file
    md_file = s_params.md_file
    try:
        root_to_kill.destroy()
    except:
        pass
    index_nr = 0
    for value in values:
        if mdp == "em":
            em_file.update(index_nr, value.get(), check_list[index_nr].get())
        elif mdp == "pr":
            pr_file.update(index_nr, value.get(), check_list[index_nr].get())
        elif mdp == "md":
            md_file.update(index_nr, value.get(), check_list[index_nr].get())
        index_nr = index_nr + 1
    save_options(em_file, pr_file, md_file, s_params)


# This function will create Simulation Steps configuration window
def steps_configure(master, restraints_button, s_params, restraintsW):
    project_name = s_params.project_name
    progress = s_params.progress
    gromacs2 = s_params.gmx_input
    if project_name != "nothing":
        root = Toplevel(master)
        root.wm_title("Simulation Steps Configuration")
        check_var1 = IntVar(root)
        check_var1.set(progress.to_do[0])
        check_var2 = IntVar(root)
        check_var2.set(progress.to_do[1])
        v1 = IntVar(root)
        v1.set(progress.x2top)
        check_var3 = IntVar(root)
        check_var3.set(progress.to_do[2])
        # Created empty variable check_var4 for genion
        check_var4 = IntVar(root)
        check_var4.set(progress.to_do[3])
        check_var5 = IntVar(root)
        check_var5.set(progress.to_do[4])
        check_var6 = IntVar(root)
        check_var6.set(progress.to_do[5])
        check_var7 = IntVar(root)
        check_var7.set(progress.to_do[6])
        check_var8 = IntVar(root)
        check_var8.set(progress.to_do[7])
        check_var9 = IntVar(root)
        check_var9.set(progress.to_do[8])
        check_var10 = IntVar(root)
        check_var10.set(progress.to_do[9])
        # Variable for Resume Simulation
        check_var11 = IntVar(root)
        check_var11.set(progress.resume)

        frame1 = Frame(root)
        frame1.pack(side=TOP)

        c1 = Checkbutton(frame1, text="Save configuration files" + steps_status_done(0, s_params), variable=check_var1,
                         command=lambda: progress.to_do_update(0, check_var1.get()))
        c1.pack(side=TOP, anchor=W)

        c2 = Checkbutton(frame1, text="Generate topology file from pdb" + steps_status_done(1, s_params),
                         variable=check_var2,
                         command=lambda: progress.to_do_update(1, check_var2.get()))
        c2.pack(side=TOP, anchor=W)

        r1 = Radiobutton(frame1, text="Use pdb2gmx tool", value=0, variable=v1,
                         command=lambda: progress.x2top_update(v1.get()))
        r1.pack(side=TOP, anchor=W)

        r2 = Radiobutton(frame1, text="Use x2top tool", value=1, variable=v1,
                         command=lambda: progress.x2top_update(v1.get()))
        r2.pack(side=TOP, anchor=W)

        c3 = Checkbutton(frame1, text="Adding Water Box (only for explicit solvent)" + steps_status_done(2, s_params),
                         variable=check_var3, command=lambda: progress.to_do_update(2, check_var3.get()))
        c3.pack(side=TOP, anchor=W)

        c4 = Checkbutton(frame1,
                         text="Adding ions and neutralize (only for explicit solvent; Optional)" + steps_status_done(3,
                                                                                                                     s_params),
                         variable=check_var4, command=lambda: progress.to_do_update(3, check_var4.get()))
        c4.pack(side=TOP, anchor=W)

        c5 = Checkbutton(frame1, text="Energy Minimization (optional)" + steps_status_done(4, s_params),
                         variable=check_var5,
                         command=lambda: progress.to_do_update(4, check_var5.get()))
        c5.pack(side=TOP, anchor=W)

        c6 = Checkbutton(frame1,
                         text="Position Restrained MD (optional, only for explicit solvent)" + steps_status_done(5,
                                                                                                                 s_params),
                         variable=check_var6, command=lambda: progress.to_do_update(5, check_var6.get()))
        c6.pack(side=TOP, anchor=W)

        c7 = Checkbutton(frame1, text="Restraints (optional)" + steps_status_done(6, s_params), variable=check_var7,
                         command=lambda: restraintsW.check(check_var7.get(), restraints_button))
        c7.pack(side=TOP, anchor=W)

        c8 = Checkbutton(frame1, text="Molecular Dynamics Simulation" + steps_status_done(7, s_params),
                         variable=check_var8,
                         command=lambda: progress.to_do_update(7, check_var8.get()))
        c8.pack(side=TOP, anchor=W)

        c9 = Checkbutton(frame1, text="Generate multimodel PDB" + steps_status_done(8, s_params), variable=check_var9,
                         command=lambda: progress.to_do_update(8, check_var9.get()))
        c9.pack(side=TOP, anchor=W)

        c10 = Checkbutton(frame1, text="Calculate vectors using ProDy (optional)" + steps_status_done(9, s_params),
                          variable=check_var10, command=lambda: progress.to_do_update(9, check_var10.get()))
        c10.pack(side=TOP, anchor=W)

        if not prody:
            check_var11.set(0)
            c10.configure(state=DISABLED)
            progress.to_do_update(9, 0)

        if gromacs2.explicit != 1:
            check_var3.set(0)
            c3.configure(state=DISABLED)
            progress.to_do_update(2, 0)
            check_var4.set(0)
            c4.configure(state=DISABLED)
            progress.to_do_update(3, 0)
            check_var6.set(0)
            c6.configure(state=DISABLED)
            progress.to_do_update(5, 0)

        l1 = Label(frame1, text="Simulation Progress:")
        l1.pack(side=TOP)

        variable_list = [check_var1, check_var2, check_var3, check_var4, check_var5, check_var6, check_var7, check_var8,
                         check_var9, check_var10, check_var11]
        progress_bar = Progressbar(frame1)
        progress_bar.pack(side=TOP)
        if check_var11.get() == 1:
            percent = steps_status_bar(check_var11.get(), s_params, variable_list)
            progress_bar.configure(value=percent)

        c11 = Checkbutton(frame1, text="Resume Simulation", variable=check_var11,
                          command=lambda: steps_click_resume(check_var11.get(), progress_bar, s_params, variable_list))
        c11.pack(side=TOP, anchor=W)

        b1 = Button(root, text="OK", command=lambda: steps_click_ok(root, s_params))
        b1.pack(side=TOP)
    elif project_name == "nothing":
        no_molecule_warning()


# This function will update status bar if checkbutton is clicked
def steps_click_resume(var, bar, s_params, variable_list=[]):
    percent = steps_status_bar(var, s_params, variable_list)
    bar.configure(value=percent)


# This function will close steps window and update number of steps to do
def steps_click_ok(root, s_params):
    root.destroy()
    progress = s_params.progress
    progress.steps = sum(progress.to_do)


# This function will show current progress on Progress Bar and operate with Steps Simulation Window for
# "Resume Simulation" button.
def steps_status_bar(var, s_params, variable_list=[]):
    progress = s_params.progress
    percent = 0.0
    if var == 1:
        to_do_nr = 0
        for step in progress.status:
            if step == 1:
                progress.to_do[to_do_nr] = 0
                progress.to_do = progress.to_do
                variable_list[to_do_nr].set(0)
            elif step == 0 and to_do_nr != 6:
                progress.to_do[to_do_nr] = 1
                progress.to_do = progress.to_do
                variable_list[to_do_nr].set(1)
            to_do_nr = to_do_nr + 1
        progress.resume = 1
    elif var == 0:
        percent = 0.0
        progress.to_do = [1, 1, 1, 1, 1, 1, 0, 1, 1, 1]
        to_do_nr = 0
        for variable in variable_list:
            if to_do_nr != 5:
                variable.set(1)
            elif to_do_nr != 5:
                variable.set(0)
            to_do_nr = to_do_nr + 1
        progress.resume = 0

    if progress.steps != 0:
        percent = ((progress.steps - sum(progress.to_do)) * 100) / progress.steps
    else:
        percent = 100

    return percent


# Gather all water options windows in one class
class WaterWindows:
    implicit_buttons = []
    explicit_buttons = []

    # Water chooser window
    def choose(self, v4_water, water_v, waterbox_button, master, s_params):
        gromacs = s_params.gmx_output
        gromacs2 = s_params.gmx_input
        root = Toplevel(master)
        root.wm_title("Water Model")

        v1 = IntVar(root)
        v1.set(gromacs2.explicit)

        v2 = IntVar(root)
        v2.set(0)

        radio_button2 = Radiobutton(root, text="Explicit Solvent Simulation", value=1, variable=v1,
                                    command=lambda: self.change_e(v1.get(), v4_water, water_v, v2, s_params))
        radio_button2.pack(side=TOP, anchor=W)

        frame1 = Frame(root, padx=10)
        frame1.pack(anchor=W)

        self.explicit_buttons = []

        for water in gromacs.water_list:
            radio_button1 = Radiobutton(frame1, text=water[1], value=water[0], variable=v4_water,
                                        command=lambda: self.change(v4_water, water_v, s_params))
            radio_button1.pack(side=TOP, anchor=W)
            self.explicit_buttons.append(radio_button1)
        self.explicit_buttons.append(waterbox_button)

        radio_button2 = Radiobutton(root, text="Implicit Solvent Simulation", value=0, variable=v1,
                                    command=lambda: self.change_e(v1.get(), v4_water, water_v, v2, s_params))
        radio_button2.pack(side=TOP, anchor=W)

        frame2 = Frame(root, padx=10)
        frame2.pack(anchor=W)
        radio_button3_1 = Radiobutton(frame2, text="Still", value=0, variable=v2,
                                      command=lambda: self.change_i(v2, s_params))
        radio_button3_1.pack(side=TOP, anchor=W)
        radio_button3_2 = Radiobutton(frame2, text="Hawkins-Cramer-Truhlar", value=1, variable=v2,
                                      command=lambda: self.change_i(v2, s_params))
        radio_button3_2.pack(side=TOP, anchor=W)
        radio_button3_3 = Radiobutton(frame2, text="Onufriev-Bashford-Case", value=2, variable=v2,
                                      command=lambda: self.change_i(v2))
        radio_button3_3.pack(side=TOP, anchor=W)

        self.implicit_buttons = [radio_button3_1, radio_button3_2, radio_button3_3]
        self.change_e(gromacs2.explicit, v4_water, water_v, v2)

        ok_button = Button(root, text="OK", command=root.destroy)
        ok_button.pack(side=TOP)

    # This function will change force field and water model when choosing Force Field in Main Window and also change
    # water model after choosing one in "waterChoose"
    def change(self, v4_water, water_v, s_params, force=False):
        gromacs = s_params.gmx_output
        gromacs2 = s_params.gmx_input
        if not force:
            force = gromacs2.force
        else:
            gromacs2.force = force
        gromacs.water_update(force)
        if gromacs2.explicit == 1:
            water_v.set(gromacs.water_list[v4_water.get() - 1][1])
        elif gromacs2.explicit == 0:
            water_v.set("Implicit Solvent")
        gromacs2.water = v4_water.get()

    # This function changes explicit to implicit and vice versa water model
    def change_e(self, value, v4_water, water_v, v2, s_params):
        progress = s_params.progress
        gromacs2 = s_params.gmx_input
        em_file = s_params.em_file
        md_file = s_params.md_file
        dynamics_dir = get_dynamics_dir()
        gromacs2.update({"explicit": value})
        if gromacs2.explicit == 1:
            for button in self.implicit_buttons:
                button.configure(state=DISABLED)
            for button in self.explicit_buttons:
                button.configure(state=ACTIVE)
            progress.to_do[2] = 1
            progress.to_do[3] = 1
            progress.to_do[5] = 1
            # em update
            if not os.path.isfile(dynamics_dir + "em.mdp"):
                parameter_nr = 0
                for parameter in em_file.options:
                    if (parameter[0] == "rlist") or (parameter[0] == ";rlist"):
                        em_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "rcoulomb") or (parameter[0] == ";rcoulomb"):
                        em_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "rvdw") or (parameter[0] == ";rvdw"):
                        em_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "implicit-solvent") or (parameter[0] == ";implicit-solvent"):
                        em_file.update(parameter_nr, "no")
                    elif (parameter[0] == "pbc") or (parameter[0] == ";pbc"):
                        em_file.update(parameter_nr, "no", 0)
                    elif (parameter[0] == "rgbradii") or (parameter[0] == ";rgbradii"):
                        em_file.update(parameter_nr, "0", 0)
                    elif (parameter[0] == "cutoff-scheme") or (parameter[0] == ";cutoff-scheme"):
                        em_file.update(parameter_nr, "Verlet")
                    elif (parameter[0] == "coulombtype") or (parameter[0] == ";coulombtype"):
                        em_file.update(parameter_nr, "PME")
                    parameter_nr = parameter_nr + 1
            # md update
            if not os.path.isfile(dynamics_dir + "md.mdp"):
                parameter_nr = 0
                for parameter in md_file.options:
                    if (parameter[0] == "nstlist") or (parameter[0] == ";nstlist"):
                        md_file.update(parameter_nr, "10")
                    elif (parameter[0] == "rlist") or (parameter[0] == ";rlist"):
                        md_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "rcoulomb") or (parameter[0] == ";rcoulomb"):
                        md_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "rvdw") or (parameter[0] == ";rvdw"):
                        md_file.update(parameter_nr, "1.0")
                    elif (parameter[0] == "Tcoupl") or (parameter[0] == ";Tcoupl"):
                        md_file.update(parameter_nr, "v-rescale")
                    elif (parameter[0] == "tau_t") or (parameter[0] == ";tau_t"):
                        md_file.update(parameter_nr, "0.1 0.1")
                    elif (parameter[0] == "tc-grps") or (parameter[0] == ";tc-grps"):
                        md_file.update(parameter_nr, "protein Non-Protein")
                    elif (parameter[0] == "ref_t") or (parameter[0] == ";ref_t"):
                        md_file.update(parameter_nr, "298 298")
                    elif (parameter[0] == "implicit-solvent") or (parameter[0] == ";implicit-solvent"):
                        md_file.update(parameter_nr, "no")
                    elif (parameter[0] == "pbc") or (parameter[0] == ";pbc"):
                        md_file.update(parameter_nr, "no", 0)
                    elif (parameter[0] == "rgbradii") or (parameter[0] == ";rgbradii"):
                        md_file.update(parameter_nr, "0", 0)
                    elif (parameter[0] == "comm_mode") or (parameter[0] == ";comm_mode"):
                        md_file.update(parameter_nr, "ANGULAR", 0)
                    elif (parameter[0] == "cutoff-scheme") or (parameter[0] == ";cutoff-scheme"):
                        md_file.update(parameter_nr, "Verlet")
                    elif (parameter[0] == "coulombtype") or (parameter[0] == ";coulombtype"):
                        md_file.update(parameter_nr, "PME")
                    parameter_nr = parameter_nr + 1
        elif gromacs2.explicit == 0:
            for button in self.implicit_buttons:
                button.configure(state=ACTIVE)
            for button in self.explicit_buttons:
                button.configure(state=DISABLED)
            progress.to_do[2] = 0
            progress.to_do[3] = 0
            progress.to_do[5] = 0
            # em update
            if not os.path.isfile(dynamics_dir + "em.mdp"):
                parameter_nr = 0
                for parameter in em_file.options:
                    if (parameter[0] == "rlist") or (parameter[0] == ";rlist"):
                        em_file.update(parameter_nr, "0")
                    elif (parameter[0] == "rcoulomb") or (parameter[0] == ";rcoulomb"):
                        em_file.update(parameter_nr, "0")
                    elif (parameter[0] == "rvdw") or (parameter[0] == ";rvdw"):
                        em_file.update(parameter_nr, "0")
                    elif (parameter[0] == "implicit-solvent") or (parameter[0] == ";implicit-solvent"):
                        em_file.update(parameter_nr, "GBSA")
                    elif (parameter[0] == "pbc") or (parameter[0] == ";pbc"):
                        em_file.update(parameter_nr, "no")
                    elif (parameter[0] == "rgbradii") or (parameter[0] == ";rgbradii"):
                        em_file.update(parameter_nr, "0")
                    elif (parameter[0] == "cutoff-scheme") or (parameter[0] == ";cutoff-scheme"):
                        em_file.update(parameter_nr, "group")
                    elif (parameter[0] == "coulombtype") or (parameter[0] == ";coulombtype"):
                        em_file.update(parameter_nr, "Cut-off")
                    parameter_nr = parameter_nr + 1
            # md update
            if not os.path.isfile(dynamics_dir + "md.mdp"):
                parameter_nr = 0
                for parameter in md_file.options:
                    if (parameter[0] == "nstlist") or (parameter[0] == ";nstlist"):
                        md_file.update(parameter_nr, "0")
                    elif (parameter[0] == "rlist") or (parameter[0] == ";rlist"):
                        md_file.update(parameter_nr, "0")
                    elif (parameter[0] == "rcoulomb") or (parameter[0] == ";rcoulomb"):
                        md_file.update(parameter_nr, "0")
                    elif (parameter[0] == "rvdw") or (parameter[0] == ";rvdw"):
                        md_file.update(parameter_nr, "0")
                    elif (parameter[0] == "Tcoupl") or (parameter[0] == ";Tcoupl"):
                        md_file.update(parameter_nr, "berendsen", 0)
                    elif (parameter[0] == "tau_t") or (parameter[0] == ";tau_t"):
                        md_file.update(parameter_nr, "0.1 0.1", 0)
                    elif (parameter[0] == "tc-grps") or (parameter[0] == ";tc-grps"):
                        md_file.update(parameter_nr, "protein Non-Protein", 0)
                    elif (parameter[0] == "ref_t") or (parameter[0] == ";ref_t"):
                        md_file.update(parameter_nr, "298 298", 0)
                    elif (parameter[0] == "implicit-solvent") or (parameter[0] == ";implicit-solvent"):
                        md_file.update(parameter_nr, "GBSA")
                    elif (parameter[0] == "pbc") or (parameter[0] == ";pbc"):
                        md_file.update(parameter_nr, "no")
                    elif (parameter[0] == "rgbradii") or (parameter[0] == ";rgbradii"):
                        md_file.update(parameter_nr, "0")
                    elif (parameter[0] == "comm_mode") or (parameter[0] == ";comm_mode"):
                        md_file.update(parameter_nr, "ANGULAR")
                    elif (parameter[0] == "cutoff-scheme") or (parameter[0] == ";cutoff-scheme"):
                        md_file.update(parameter_nr, "group")
                    elif (parameter[0] == "coulombtype") or (parameter[0] == ";coulombtype"):
                        md_file.update(parameter_nr, "Cut-off")
                    parameter_nr = parameter_nr + 1
            self.change_i(v2, s_params)
            # in implicit solvent watermodel must be set to "None"
            v4_water.set(len(self.explicit_buttons) - 1)

        self.change(v4_water, water_v, s_params)

    # This function changes implicit water model
    @staticmethod
    def change_i(int_variable, s_params):
        em_file = s_params.em_file
        md_file = s_params.md_file
        dynamics_dir = get_dynamics_dir()
        if int_variable.get() == 0:
            if not os.path.isfile(dynamics_dir + "em.mdp"):
                parameter_nr = 0
                for parameter in em_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        em_file.update(parameter_nr, "Still")
                    parameter_nr = parameter_nr + 1
            if not os.path.isfile(dynamics_dir + "md.mdp"):
                parameter_nr = 0
                for parameter in md_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        md_file.update(parameter_nr, "Still")
                    parameter_nr = parameter_nr + 1
        elif int_variable.get() == 1:
            if not os.path.isfile(dynamics_dir + "em.mdp"):
                parameter_nr = 0
                for parameter in em_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        em_file.update(parameter_nr, "HCT")
                    parameter_nr = parameter_nr + 1
            if not os.path.isfile(dynamics_dir + "md.mdp"):
                parameter_nr = 0
                for parameter in md_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        md_file.update(parameter_nr, "HCT")
                    parameter_nr = parameter_nr + 1
        elif int_variable.get() == 2:
            if not os.path.isfile(dynamics_dir + "em.mdp"):
                parameter_nr = 0
                for parameter in em_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        em_file.update(parameter_nr, "OBC")
                    parameter_nr = parameter_nr + 1
            if not os.path.isfile(dynamics_dir + "md.mdp"):
                parameter_nr = 0
                for parameter in md_file.options:
                    if (parameter[0] == "gb-algorithm") or (parameter[0] == ";gb-algorithm"):
                        md_file.update(parameter_nr, "OBC")
                    parameter_nr = parameter_nr + 1

    # Water box configuration window
    @staticmethod
    def box(master, s_params):
        gromacs2 = s_params.gmx_input
        root = Toplevel(master)
        root.wm_title("Water Box Options")
        root.wm_geometry("300x200")
        v = StringVar(root)
        v.set(gromacs2.box_type)
        w = Label(root, text="Box type")
        w.pack()
        radio_button = Radiobutton(root, text="triclinic", value="triclinic", variable=v,
                                   command=lambda: gromacs2.update({"box_type": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="cubic", value="cubic", variable=v,
                                   command=lambda: gromacs2.update({"box_type": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="dodecahedron", value="dodecahedron", variable=v,
                                   command=lambda: gromacs2.update({"box_type": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="octahedron", value="octahedron", variable=v,
                                   command=lambda: gromacs2.update({"box_type": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        w1 = Label(root, text="Distance")
        w1.pack()
        distance = Entry(root)
        distance.pack(side=TOP)
        distance.insert(0, gromacs2.box_distance)
        w2 = Label(root, text="Density [g/L]")
        w2.pack()
        density = Entry(root)
        density.pack(side=TOP)
        density.insert(0, gromacs2.box_density)
        ok_button = Button(root, text="OK", command=lambda: gromacs2.update(
            {"box_distance": distance.get(), "box_density": density.get()}, root))
        ok_button.pack(side=TOP)

    # Hydrogen configuration (for bigger time steps)
    @staticmethod
    def box2(master, s_params):
        gromacs2 = s_params.gmx_input
        root = Toplevel(master)
        root.wm_title("Hydrogen options (for Pdb2gmx)")
        root.wm_geometry("300x200")
        v1 = StringVar(root)
        v1.set(gromacs2.hydro)
        w = Label(root, text="Hydrogen type (for pdb2gmx only)")
        w.pack()
        radio_button = Radiobutton(root, text="Normal Hydrogen", value="noheavyh", variable=v1,
                                   command=lambda: gromacs2.update({"hydro": v1.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="Deuterium", value="deuterate", variable=v1,
                                   command=lambda: gromacs2.update({"hydro": v1.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="Heavy Hydrogen (4amu) ", value="heavyh", variable=v1,
                                   command=lambda: gromacs2.update({"hydro": v1.get()}))
        radio_button.pack(side=TOP, anchor=W)
        ok_button = Button(root, text="OK", command=root.destroy)
        ok_button.pack(side=TOP)


# Options for the genion class all the options
class GenionWindow:

    # Genion box configuration window
    def window(self, master, s_params):
        gromacs2 = s_params.gmx_input
        root = Toplevel(master)
        root.wm_title("GENION options")
        root.wm_geometry("300x350")
        v = StringVar(root)
        v.set(gromacs2.neutrality)
        w = Label(root, text="Parameters for genion")
        w.pack()
        radio_button = Radiobutton(root, text="Neutralize System", value="neutral", variable=v,
                                   command=lambda: gromacs2.update({"neutrality": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        radio_button = Radiobutton(root, text="Do not Neutralize", value="noneutral", variable=v,
                                   command=lambda: gromacs2.update({"neutrality": v.get()}))
        radio_button.pack(side=TOP, anchor=W)
        w1 = Label(root, text="Salt Concentration")
        w1.pack()
        salt = Entry(root)
        salt.pack(side=TOP)
        salt.insert(0, gromacs2.salt_conc)
        w2 = Label(root, text="Positive Ion")
        w2.pack()
        posit = Entry(root)
        posit.pack(side=TOP)
        posit.insert(0, gromacs2.positive_ion)
        w3 = Label(root, text="Negative Ion")
        w3.pack()
        negat = Entry(root)
        negat.pack(side=TOP)
        negat.insert(0, gromacs2.negative_ion)
        ok_button = Button(root, text="OK", command=lambda: gromacs2.update(
            {"salt_conc": salt.get(), "positive_ion": posit.get(), "negative_ion": negat.get()}, root))
        ok_button.pack(side=TOP)

# This is the window to setup ProDy options
# def vectors_window(master, s_params):
#    project_name = s_params.project_name
#    if project_name != "nothing":
#        root = Toplevel(master)
#        root.wm_title("Vectors Configuration")
#
#        frame1 = Frame(root)
#        frame1.pack()
#
#        v1 = IntVar(root)
#        v1.set(calculation_type)
#        v2 = IntVar(root)
#        v2.set(contact_map)
#
#        radio_button0 = Radiobutton(frame1, text="Anisotropic network model", value=0, variable=v1,
#                                        command=lambda: block_contact(0, c1, v2))
#        radio_button0.pack()
#        radio_button1 = Radiobutton(frame1, text="Principal component analysis", value=1, variable=v1,
#                                        command=lambda: block_contact(1, c1, v2))
#        radio_button1.pack()
#        radio_button2 = Radiobutton(frame1, text="Gaussian network model (experimental)", value=2, variable=v1,
#                                        command=lambda: block_contact(0, c1, v2))
#        radio_button2.pack()
#
#        c1 = Checkbutton(frame1, text="Show Contact Map", variable=v2)
#        c1.pack()
#        if block_contact_map == 1:
#            c1.configure(state=DISABLED)
#
#        ok_button = Button(frame1, text="OK", command=lambda: options_change(v1, v2, root))
#        ok_button.pack(side=TOP)
#
#    elif project_name == "nothing":
#        no_molecule_warning()
