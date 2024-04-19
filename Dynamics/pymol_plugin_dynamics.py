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
# - Natalia Floriańska

from __future__ import print_function

# --Import libraries--
# Import native python libraries
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
from tkinter import *
from tkinter import messagebox as tkMessageBox
from tkinter import filedialog as tkFileDialog
from tkinter.ttk import Progressbar, Scrollbar
# Import libraries from PyMOL specific work.
from pymol import cmd, parsing, plugins, CmdException

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

# Import Dynamics libraries
from version import VERSION
import SimulationParameters
import MdpConfig
import CalculationWindow
import InterpretationWindow
import WaterWindows
import RestraintsWindow
import GenionWindow

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
calculationW = CalculationWindow.CalculationWindow()


# This function will initialize all plugin stuffs
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
    simulation_parameters = SimulationParameters.SimulationParameters()

    supported_gmx_versions = ["2016", "2018"]
    if not len(simulation_parameters.gmx_output.gmx_exe):
        print("GROMACS 2016 or newer not detected.")
        status = ["fail",
                  "GROMACS not detected. Please install and setup GROMACS 2016 or newer correctly for your platform."
                  " Check '~/.dynamics/test_gromacs.txt' for more details. Don't forget to add GROMACS bin directory"
                  " to your PATH"]
    elif simulation_parameters.gmx_output.version[0:4] not in supported_gmx_versions:
        print("Warning. Unsupported GROMACS Version")
    if not travis_ci:
        create_gui(gui_library, status, simulation_parameters, parent)
    return status, simulation_parameters


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
        em_file = MdpConfig.MdpConfig("em.mdp", em_file_config, 1)
    else:
        em_file = MdpConfig.MdpConfig("em.mdp", EM_INIT_CONFIG, 0)
    if os.path.isfile(dynamics_dir + "pr.mdp"):
        shutil.copy(dynamics_dir + "pr.mdp", project_dir + "pr.mdp")
        print("Found pr.mdp file. Using it instead of local configuration.")
    elif os.path.isfile(project_dir + "pr.mdp"):
        pr_file_config = open(project_dir + "pr.mdp", "r").read()
        pr_file = MdpConfig.MdpConfig("pr.mdp", pr_file_config, 1)
    else:
        pr_file = MdpConfig.MdpConfig("pr.mdp", PR_INIT_CONFIG, 0)
    if os.path.isfile(dynamics_dir + "md.mdp"):
        shutil.copy(dynamics_dir + "md.mdp", project_dir + "md.mdp")
        print("Found md.mdp file. Using it instead of local configuration.")
    elif os.path.isfile(project_dir + "md.mdp"):
        md_file_config = open(project_dir + "md.mdp", "r").read()
        md_file = MdpConfig.MdpConfig("md.mdp", md_file_config, 1)
    else:
        md_file = MdpConfig.MdpConfig("md.mdp", MD_INIT_CONFIG, 0)
    #    save_options()
    try:
        if project_name in cmd.get_names("objects"):  # PyMOL API
            cmd.save(project_dir + project_name + ".pdb", project_name)  # PyMOL API
            print("cmd saved")
    except (AttributeError, TypeError) as e:
        pass
    return em_file, pr_file, md_file


def get_dynamics_dir():
    home_dir = os.path.expanduser('~')
    gmx_home_dir_path = os.path.abspath(home_dir)
    dynamics_dir = os.path.join(gmx_home_dir_path, '.dynamics', '')
    return dynamics_dir


def get_project_dirs(project_name="nothing"):
    dynamics_dir = get_dynamics_dir()
    project_dir = os.path.join(dynamics_dir, project_name, '')
    return project_dir


def check_output_subprocess(command, stdin_file_path=None):
    stdin_file = None
    stdin_msg = "None"
    if stdin_file_path:
        stdin_file = open(stdin_file_path, "r")
        stdin_msg = stdin_file_path
    print("Running command: " + command + "; STDIN: " + stdin_msg)
    output = subprocess.check_output(command, shell=True, stdin=stdin_file).decode(sys.stdout.encoding)
    if stdin_file_path:
        stdin_file.close()

    return output


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


# Read lines of text file and standardize new line char. Returns list.
def read_lines_standardize(text_file_path):
    text_lines = []

    ifs = open(text_file_path, "r")
    for line in iter(ifs.readline, ''):
        new_line = re.sub("(\r\n)|(\r)", "\n", line)
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
    if( len(waters_info)<=0):
        waters_info = get_alternate_water_info(gmx_output_lines)
    waters_info2 = []
    number = 1
    for water in waters_info:
        waters_info2.append([number, water[:-1]])
        number = number + 1

    return waters_info2
#In newwer version of Gromacs water is under Section Select Water Model
def get_alternate_water_info(gmx_output_lines):
    wm_start_line = 0
    while gmx_output_lines[wm_start_line] != "Select the Water Model:\n":
            wm_start_line = wm_start_line + 1
    wm_start_line = wm_start_line + 2
    wm_end_line = wm_start_line
    while gmx_output_lines[wm_end_line] != "\n":
            wm_end_line = wm_end_line + 1
    return gmx_output_lines[wm_start_line:wm_end_line]


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
    print("Status after update")
    print(status)


# This function will start real workflow of the plugin, once everything is set
def dynamics(s_params):
    print("Starting PyMOL plugin 'dynamics' ver. {}".format(plugin_ver))    
    status = ["ok", ""]
    project_name = s_params.project_name
    project_dir = get_project_dirs(project_name)
    progress = s_params.progress
    gromacs2 = s_params.gmx_input
    vectors_prody = s_params.vectors_prody
    print("Vectors prody")
    print(vectors_prody)
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
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    elif status[0] == "ok" and stop == 0 and progress.to_do[1] == 1 and progress.x2top == 1:
        status = gromacs2.x2top(s_params)        
        if status[0] == "ok":
            progress.status[1] = 1
            progress.to_do[1] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    # Adding water box
    if status[0] == "ok" and stop == 0 and progress.to_do[2] == 1:
        status = gromacs2.waterbox(s_params)        
        if status[0] == "ok":
            progress.status[2] = 1
            progress.to_do[2] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    # Adding ions
    if status[0] == "ok" and stop == 0 and progress.to_do[3] == 1:
        status = gromacs2.saltadd(s_params)        
        if status[0] == "ok":
            progress.status[3] = 1
            progress.to_do[3] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    # EM
    if status[0] == "ok" and stop == 0 and progress.to_do[4] == 1:
        status = gromacs2.em(s_params)        
        if status[0] == "ok":
            progress.status[4] = 1
            progress.to_do[4] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)
    elif status[0] == "ok" and stop == 0 and progress.to_do[4] == 0 and progress.status[4] == 0:
        shutil.copy(project_name + "_b4em.gro", project_name + "_b4pr.gro")

    # PR
    if status[0] == "ok" and stop == 0 and progress.to_do[5] == 1:
        status = gromacs2.pr(s_params)        
        if status[0] == "ok":
            progress.status[5] = 1
            progress.to_do[5] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)
    elif status[0] == "ok" and stop == 0 and progress.to_do[5] == 0 and progress.status[5] == 0:
        shutil.copy(project_name + "_b4pr.gro", project_name + "_b4md.gro")

    # Restraints
    if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
        status = gromacs2.restraints(project_name)        
        if status[0] == "ok":
            progress.status[6] = 1
            progress.to_do[6] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    # MD
    if status[0] == "ok" and stop == 0 and progress.to_do[7] == 1:
        status = gromacs2.md(s_params)
        if status[0] == "ok":
            progress.status[7] = 1
            progress.to_do[7] = 0
            calculationW.bar_update(s_params, status)
            save_options(s_params)

    # Trjconv
    if status[0] == "ok" and stop == 0 and progress.to_do[8] == 1:
        status = gromacs2.trjconv(s_params)
        show_multipdb(s_params)
        progress.status[8] = 1
        progress.to_do[8] = 0
        calculationW.bar_update(s_params, status)
        save_options(s_params)

    # Calculating vectors
    if status[0] == "ok" and stop == 0 and progress.to_do[9] == 1 and prody:
        vectors_prody.prody(project_name)
        vectors_prody.nmd_format(project_name)
        vectors_prody.show_vectors()
        progress.status[9] = 1
        progress.to_do[9] = 0
        calculationW.bar_update(s_params, status)
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
    print("Load Options 1")
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

    if options[0][1:4] == "2.2" or options[0][1:4] == "3.0":
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
    waterW = WaterWindows.WaterWindows()
    restraintsW = RestraintsWindow.RestraintsWindow()
    genionW = GenionWindow.GenionWindow()
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
                                                                              v4_water, water_v, check1_button, s_params))
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






# Show interpretation window...
def show_interpretation_window(parent, s_params):
    InterpretationWindow.InterpretationWindow(parent, s_params)


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
        s_params.project_name = project_name
        s_params.create_cfg_files()
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
    if file:
        load_file(file.name, s_params)
        v1_name.set(project_name)
        v2_group.set(gromacs.group_list[gromacs2.group][0])
        v3_force.set(gromacs.force_list[gromacs2.force - 1][0])
        v4_water.set(gromacs.water_list[gromacs2.water - 1][0])
        water_v.set(gromacs.water_list[v4_water.get() - 1][1])
        radio_button1 = Radiobutton(frame1_1a, text=project_name, value=project_name, variable=v1_name,
                                    command=lambda: set_variables(v1_name.get(), v2_group, v3_force, v4_water, water_v,
                                                                  config_button_restraints, s_params))
        radio_button1.pack(side=TOP, anchor=W)
    root.destroy()


# This function sets variables after choosing new molecule
def set_variables(name, v2_group, v3_force, v4_water, water_v, config_button_restraints, s_params):
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
        s_params.project_name = project_name
        s_params.create_cfg_files()
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


# This function will update MDP class objects after closing "mdp_configure" window
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
    save_options(s_params)


# This function will create Simulation Steps configuration window
def steps_configure(master, restraints_button, s_params, restraints_w):
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
                         command=lambda: restraints_w.check(check_var7.get(), restraints_button))
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


# This function will close steps window and update number of steps to take
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
        # print("progress.steps: "+str(progress.steps))
        # print("progress.to_do: "+str(progress.to_do))
        # percent = ((progress.steps - sum(progress.to_do)) * 100) / progress.steps
        # percent = sum(progress.to_do) / progress.steps * 100
        percent = sum(progress.status)/sum(progress.to_do) * 100 
    else:
        percent = 100

    return percent


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
