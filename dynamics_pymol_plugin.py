#!/usr/bin/env python2
#-*- coding: utf-8 -*-

##This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
##Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
##Contributors:
##- Tomasz Makarewicz (tomaszm@biotech.ug.edu.pl)

##--Import libraries--
##Import nativ python libraries
import subprocess, time, os, shutil, thread, pickle
##Import libraries for tk graphic interface
from Tkinter import *
from Tix import *
import tkSimpleDialog, tkMessageBox, tkFileDialog
##Import libraries from PyMOL specific work. Detect if running as a plugin. If true set plugin variable to 1.
try:
	from pymol import cmd, stored
	cmd.get_version()
	plugin = 1
except:
	plugin = 0
##Globals variables
help_name = ["-h", "h", "-help", "help"]
clean_name = ["-c", "c", "clean", "-clean"]
plugin_ver = " 1.1.0pre"

stop = 0
status = ["ok", ""]
error = ""

em_init_config = """cpp = /usr/bin/cpp
define = -DFLEX_SPC
constraints = none
integrator = steep
nsteps = 100
nstlist = 10
ns_type = grid
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
emtol = 1000.0
emstep = 0.01"""

pr_init_config = """cpp = /usr/bin/cpp
define = -DPOSRES
constraints = all-bonds
integrator = md
dt = 0.002
nsteps = 500
nstcomm = 1
nstxout = 10
nstvout = 1000
nstfout = 0
nstlog = 10
nstenergy = 10
nstlist = 10
ns_type = grid
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
Tcoupl = berendsen
tau_t = 0.1 0.1
tc-grps = protein Non-Protein
ref_t = 300 300
Pcoupl = no
tau_p = 0.5
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = yes
gen_temp = 300.0
gen_seed = 173529"""

md_init_config = """cpp = /usr/bin/cpp
;define = -DPOSRES
constraints = all-bonds
integrator = md
dt = 0.002
nsteps = 5000
nstcomm = 1
nstxout = 50
nstvout = 0
nstfout = 0
nstlist = 10
ns_type = grid
rlist = 1.0
rcoulomb = 1.0
rvdw = 1.0
Tcoupl = berendsen
tau_t = 0.1 0.1
tc-grps = protein Non-Protein
ref_t = 300 300
Pcoupl = no
tau_p = 0.5
compressibility = 4.5e-5
ref_p = 1.0
gen_vel = yes
gen_temp = 300.0
gen_seed = 173529"""

project_name = 'nothing'
dir_path_dynamics = os.getenv("HOME")+'/.dynamics/'
dir_path_project = dir_path_dynamics+project_name + '/'
##Clean "nothing" temporary directory if present.
try:
	shutil.rmtree(dir_path_project)
except:
	pass
##Creating new directories
if os.path.isdir(dir_path_project) == False:
	os.makedirs(dir_path_project)

##This class is responsible for interface to GROMACS. It will read all important data from GROMACS tools.
class Gromacs_output:
	
	version = "GROMACS not found"
	status = ["ok", ""]
	gromacs_path = ""
	force_list = []
	water_list = []
	group_list = []
	restraints = []
	
	def __init__(self):
		for_water = "1"
		status = ["ok", ""]
		gromacs_path = ""
		print "Testing GROMACS installation and version"
		subprocess.call("echo -e '"+for_water+"\n1' | pdb2gmx &> "+dir_path_dynamics+"test_gromacs.txt", executable="/bin/bash", shell=True)
		test_gromacs = open(dir_path_dynamics+"test_gromacs.txt","r")
		lista_gromacs = test_gromacs.readlines()
		#This will require correction in case GROMACS header will change
		if lista_gromacs[0] == "                         :-)  G  R  O  M  A  C  S  (-:\n":
			version = lista_gromacs[4].split("  ")
			gromacs_version = "GROMACS " + version[15]
			print "Found " + gromacs_version
		else:
			print "GROMACS not detected."
			import platform
			print "System: "+platform.platform()
			if platform.system() == "Linux":
				if platform.machine() == "i686":
					gromacs_path = "export PATH="+dir_path_dynamics+'gromacs-4.5.5-linux-32/bin:"${PATH}"; export LD_LIBRARY_PATH='+dir_path_dynamics+"/gromacs-4.5.5-linux-32/lib; "
					gromacs_version = "GROMACS 4.5.5"
					if os.path.isdir(dir_path_dynamics+"gromacs-4.5.5-linux-32/") == False:
						import urllib, tarfile
						print "Downloading GROMACS 4.5.5 for your platform"
						urllib.urlretrieve("http://ubuntuone.com/3MhyjK4mfavbCqdAzCbjA9", dir_path_dynamics+"gromacs-4.5.5-linux-32.tar.bz2")
						tarfile.open(dir_path1+"gromacs-4.5.5-linux-32.tar.bz2").extractall(dir_path_dynamics+"gromacs-4.5.5-linux-32/")
				elif platform.machine() == "x86_64":
					gromacs_path = "export PATH="+dir_path_dynamics+'gromacs-4.5.5-linux-64/bin:"${PATH}"; export LD_LIBRARY_PATH='+dir_path_dynamics+"/gromacs-4.5.5-linux-64/lib; "
					gromacs_version = "GROMACS 4.5.5"
					if os.path.isdir(dir_path_dynamics+"gromacs-4.5.5-linux-64/") == False:
						import urllib, tarfile
						print "Downloading GROMACS 4.5.5 for your platform"
						urllib.urlretrieve("http://ubuntuone.com/27vlaNtNV6mDHs7qu7qXYv", dir_path_dynamics+"gromacs-4.5.5-linux-64.tar.bz2")
						tarfile.open(dir_path_dynamics+"gromacs-4.5.5-linux-64.tar.bz2").extractall(dir_path_dynamics+"gromacs-4.5.5-linux-64/")
				else:
					print "Please install and setup correctly GROMACS for your platform. Aborting."
					status = ["fail", "Require GROMACS installation"]
			else:
				print "Please install and setup correctly GROMACS for your platform. Aborting."
				status = ["fail", "Require GROMACS installation"]
		self.status = status
		self.path = gromacs_path
		if status[0] == "ok":
			self.version = gromacs_version

		subprocess.call(self.path+"echo -e '"+for_water+"\n1' | pdb2gmx &> "+dir_path_dynamics+"test_gromacs.txt", executable="/bin/bash", shell=True)
		test_gromacs = open(dir_path_dynamics+"test_gromacs.txt","r")
		lista_gromacs = test_gromacs.readlines()
		print "Reading available force forces"	
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

		print "Reading available water models"
		water_start_line = 0
		while lista_gromacs[water_start_line][0:7] != "Opening":
			water_start_line = water_start_line + 1
		water_start_line = water_start_line + 1
		water_end_line = water_start_line + 1
		while (lista_gromacs[water_end_line][0:7] != "Opening") and (lista_gromacs[water_end_line][0] != "\n"):
			water_end_line = water_end_line + 1
		water_list = lista_gromacs[water_start_line:water_end_line]
		water_list2 = []
		number = 1
		for water in water_list:
			water_list2.append([number, water[:-1]])
			number = number + 1
	
		print "Reading available groups"
		#Temporary solution. Need real reading from GROMACS.
		group_list = ["0: System", "1: Protein", "2: Protein-H", "3: C-alpha", "4: Backbone", "5: MainChain", "6: MainChain+Cb", "7: MainChain+H", "8: SideChain", "9: SideChain-H", "10: Prot-Masses", "11: non-Protein", "12: Water", "13: SOL", "14: non-Water"]
		group_list2 = []
		number = 0
		for group in group_list:
			group_list2.append([number, group])
			number = number + 1
		
		self.force_list = force_list2
		self.water_list = water_list2
		self.group_list = group_list2
	
	##This function will update water list if force field is changed.
	def water_update(self, force_number):
		print "Updateing available water models"
		subprocess.call(self.path+"echo -e '"+str(force_number)+"\n1' | pdb2gmx &> "+dir_path_dynamics+"test_gromacs.txt", executable="/bin/bash", shell=True)
		test_gromacs = open(dir_path_dynamics+"test_gromacs.txt","r")
		lista_gromacs = test_gromacs.readlines()

		water_start_line = 0
		while lista_gromacs[water_start_line][0:7] != "Opening":
			water_start_line = water_start_line + 1
		water_start_line = water_start_line + 1
		water_end_line = water_start_line + 1
		while (lista_gromacs[water_end_line][0:7] != "Opening") and (lista_gromacs[water_end_line][0] != "\n"):
			water_end_line = water_end_line + 1
		water_list = lista_gromacs[water_start_line:water_end_line]
		water_list2 = []
		number = 1
		for water in water_list:
			water_list2.append([number, water[:-1]])
			number = number + 1
			
		self.water_list = water_list2
		save_options()
		return water_list2
	
	##This function will read atoms group for restraints for current molecule.	
	def index(self):
		self.restraints = []
		os.chdir(dir_path_project)
		subprocess.call(self.path+"echo q | make_ndx -f "+project_name+".pdb -o index.ndx &> restraints.log", executable="/bin/bash", shell=True)	
		index = open("index.ndx","r")
		index_list = index.readlines()
		index_position = 0
		atoms = ""
		for line in index_list:
			if line[0] == "[":
				self.restraints.append([])
				self.restraints[index_position].append(line)
				if index_position != 0:
					self.restraints[index_position-1].append(atoms)
				index_position = index_position + 1
				atoms = ""
			else:
				atoms = atoms + line
		self.restraints[index_position-1].append(atoms)

gromacs = Gromacs_output()

##This class is responsible for performing molecular dynamics simulation with GROMACS tools.
class Gromacs_input:
	
	status = ["ok", "Waiting to start"]
	force = 1
	water = 1
	group = 1
	box_type = "triclinic"
	box_distance = "0.5"
	box_density = "1000"
	restraints_nr = 1

	##This function will change given variabless stored by the class (needed for lambda statements)
	def update(self, group, box_type, box_distance, box_density, root=""):
		#Close mother window if present
		try:
			root.destroy()
		except:
			pass
		
		self.group = group
		self.box_type = box_type
		self.box_distance = box_distance
		self.box_density = box_density
		save_options()
		print "gromacs update"
	
	##This function will create initial topology and triectory using pdb file and choosen force force
	def pdb2top(self, file_path, force):
		status = ["ok", "Calculating topology using Force forces"]
		self.status = status
		try:
			os.remove(project_name+".gro")
			os.remove(project_name+".top")
		except:
			pass
		print "Calculating topology using Force forces"
		Pdb2gmx = subprocess.call(gromacs.path+"echo -e '"+force+"' | pdb2gmx -f "+project_name+".pdb -o "+project_name+".gro -p "+project_name+".top &> log.txt",
		executable="/bin/bash", shell=True)

		if os.path.isfile(file_path+".gro") == True:
			status = ["ok", ""]
		else:
			status = ["fail", "Warning. Trying to ignore unnecessary hydrogen atoms."]
			print status[1]
			Pdb2gmx = subprocess.call(gromacs.path+"echo -e '"+force+"' | pdb2gmx -ignh -f "+project_name+".pdb -o "+project_name+".gro -p "+project_name+".top &> log.txt",
			executable="/bin/bash", shell=True)

		if os.path.isfile(file_path+".gro") == True:
			status = ["ok", "Calculated topology using Force forces"]
		else:
			status = ["fail", "Force force unable to create topology file"]
		
		self.status = status
	
	##This function will create and add waterbox.
	def waterbox(self, file_path):
		status = ["ok", "Adding Water Box"]
		self.status = status
		box_type = "-bt "+self.box_type+" "
		distance = "-d "+self.box_distance+" "
		density = "-density "+self.box_density
		
		try:
			os.remove(project_name+"1.gro")
		except:
			pass
		
		print "Generating water_box"
		Editconf = subprocess.call(gromacs.path+"editconf -f "+project_name+".gro -o "+project_name+"1.gro "+box_type+distance+density+" &>> log.txt",
		executable="/bin/bash", shell=True)

		print "Adding Water Box"
		Genbox = subprocess.call(gromacs.path+"genbox -cp "+project_name+"1.gro -cs -o "+project_name+"_b4em.gro -p "+project_name+".top &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"1.gro") == True:
			status = ["ok", "Added Water Box"]
		else:
			status = ["fail", "Unable to add waterbox"]
		self.status = status
	
	##This function will perform energy minimalization	
	def em(self, file_path):
		status = ["ok", "Energy Minimalization"]
		self.status = status
		
		try:
			os.remove(project_name+"_em.tpr")
			os.remove(project_name+"_em.trr")
		except:
			pass

		print "Energy Minimalization"		
		Grompp = subprocess.call(gromacs.path+"grompp -f em -c "+project_name+"_b4em -p "+project_name+" -o "+project_name+"_em &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+project_name+"_em -o "+project_name+"_em -c "+project_name+"_b4pr -v &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"_em.tpr") == True:
			status = ["ok", "Energy Minimalized"]
		else:
			status = ["fail", "Unable to perform Energy Minimalization"]
		self.status = status
	
	##This function will perform position restrained MD
	def pr(self, file_path):
		status = ["ok", "Position Restrained MD"]
		self.status = status
		
		try:
			os.remove(project_name+"_pr.tpr")
			os.remove(project_name+"_pr.trr")
		except:
			pass
		
		print "Position Restrained MD"
		Grompp = subprocess.call(gromacs.path+"grompp -f pr -c "+project_name+"_b4pr -r "+project_name+"_b4pr -p "+project_name+" -o "+project_name+"_pr &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+project_name+"_pr -o "+project_name+"_pr -c "+project_name+"_b4md -v &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"_pr.tpr") == True:
			status = ["ok", "Position Restrained MD finished"]
		else:
			status = ["fail", "Unable to perform Position Restrained"]
		self.status = status
	
	##This function will create posre.itp file for molecular dynamics simulation with choosen atoms if restraints were selected
	def restraints(self):
		status = ["ok", "Adding Restraints"]
		self.status = status
		
		try:
			os.remove("posre_2.itp")
		except:
			pass
			
		print "Adding Restraints"
		Genrestr = subprocess.call("echo 0 | genrestr -f "+project_name+".pdb -o posre_2.itp -n index_dynamics.ndx &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile("posre_2.itp") == True:
			status = ["ok", "Added Restraints"]
			os.remove("posre.itp")
			shutil.copy("posre_2.itp", "posre.itp")
		else:
			status = ["fail", "Unable to create restraints file"]
		self.status = status
	
	##This function will perform position final molecular dynamics simulation
	def md(self, file_path):
		status = ["ok", "Molecular Dynamics Simulation"]
		self.status = status
		
		try:
			os.remove(project_name+"_md.tpr")
			os.remove(project_name+"_md.trr")
		except:
			pass
		
		print "Molecular Dynamics Simulation"
		Grompp = subprocess.call(gromacs.path+"grompp -f md -c "+project_name+"_b4md  -p "+project_name+" -o "+project_name+"_md &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+project_name+"_md -o "+project_name+"_md -c "+project_name+"_after_md -v &>> log.txt", executable="/bin/bash", shell=True)
	
		if os.path.isfile(file_path+"_md.tpr") == True:
			status = ["ok", "Molecular Dynamics Simulation finished"]
		else:
			status = ["fail", "Unable to perform Molecular Dynamics Simulation"]
		self.status = status
	
	##This function will convert final results to multimodel pdb file
	def trjconv(self, file_path, group):
		status = ["ok", "Creating Multimodel PDB"]
		self.status = status
		
		try:
			os.remove(project_name+"_multimodel.pdb")
		except:
			pass
		if os.path.isfile(project_name+"_multimodel.pdb") == True:
			os.remove(project_name+"_multimodel.pdb")
		
		print "Creating Multimodel PDB"
		Trjconv = subprocess.call(gromacs.path+"echo "+group+" | trjconv -f "+project_name+"_md.trr -s "+project_name+"_md.tpr -app -o "+project_name+"_multimodel.pdb &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"_multimodel.pdb") == True:
			status = ["ok", "Finished!"]
		else:
			status = ["fail", "Unable to generate multimodel PDB file"]
		self.status = status

gromacs2 = Gromacs_input()

##This class create and maintain abstraction mdp file representatives. em.mdp, pr.mdp, md.mdp
class Mdp_config:
	
	config = ""
	options = [[]]
	file_name = ""
	
	def __init__(self, file_name, init_config, external_file=0):
		if external_file == 0:
			self.config = """title = """+project_name+"_"+file_name+"\n"+init_config
		elif external_file == 1:
			self.config = init_config
		self.file_name = file_name
		list1 = self.config.split("\n")
		list2 = []
		for line in list1:
			list2.append(line.split(" = "))
		self.options = list2			
	
	def update(self, values=""):
		if values != "":
			index_nr = 0
			for option in self.options:
				option[1] = values[index_nr]
				index_nr = index_nr + 1
		config = ""
		for option in self.options:
			config = config + option[0]+" = "+option[1]+"\n"
		self.config = config
		save_options()
	
	def save_file(self):
		mdp = open(dir_path_project+self.file_name, "w")
		mdp.write(self.config)
		mdp.close() 

##Status and to_do maintaining class
class Progress_status:
	
	status= []
	to_do = []
	from_begining = 1
	
	def __init__(self, dir_path):
		if os.path.isfile(dir_path+"/status.pickle"):
			pickle_file = file(dir_path+"/status.pickle")
			self.status = pickle.load(pickle_file)
		else:
			self.status = [0,0,0,0,0,0,0]
		self.pickle_file = dir_path+"/status.pickle"
		self.to_do = [1,1,1,1,1,1,1]
	
	def to_do_update(self, position, value):
		if type(position) == type(1):
			self.to_do[position] = value
	
	def to_do_status(self):
		to_do = []
		for work in self.status:
			if work == 0:
				to_do.append(1)
			elif work == 1:
				to_do.append(0)
		self.to_do = to_do

##init function - puts plugin into menu and starts 'dynamicsDialog' after clicking.
def __init__(self):
	self.menuBar.addmenuitem("Plugin", "command", "dynamics"+plugin_ver, label = "dynamics"+plugin_ver,
	command = rootWindow)

##--Graphic Interface--
##Root menu window
def rootWindow():

	root = Tk()
	root.wm_title("Dynamics"+plugin_ver)
	
	##Creating objects
	global calculationW, restraintsW
	calculationW = CalculationWindow()
	restraintsW = RestraintsWindow()
	
	##Detect list of PyMOL loaded PDB files if no files than list "nothing"
	if plugin == 1:
		allNames = cmd.get_names("objects") #PyMOL API
		allNames1 = []
		for name in allNames:
			name1 = name.split("_")
			if name1[-1] == "multimodel" or name1[-1] == "(sele)":
				pass
			else:
				allNames1.append(name)
		allNames = allNames1
		
	elif plugin == 0:
		allNames = ["nothing"]
	if allNames == []:
		allNames = ["nothing"]
	
	##TkInter variables
	global dir_path_project, project_name
	project_name = allNames[0]
	dir_path_project = dir_path_dynamics + project_name + '/'
	
	v1_name = StringVar(root)
	v1_name.set(project_name)
	
	groupNr = gromacs.group_list[1][0]
	v2_group = IntVar(root)
	v2_group.set(groupNr)
	
	forceNr = gromacs.force_list[0][0]
	v3_force = IntVar(root)
	v3_force.set(forceNr)
	
	waterNr = gromacs.water_list[0][0]
	v4_water = IntVar(root)
	v4_water.set(waterNr)

	water_v = StringVar(root)
	water_v.set(gromacs.water_list[v4_water.get()-1][1])
	
	check_var = IntVar(root)
	
	check_var2 = IntVar(root)
	check_var2.set(restraintsW.on_off)
	
	#Initial configuration
	set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var)
	
	##Start drawing interface
	frame0 = Frame(root)
	frame0.pack(side=TOP)
	
	w_version = Label(frame0, text=gromacs.version)
	w_version.pack(side=TOP)
	
	frame1 = Frame(root)
	frame1.pack(side=TOP)
	
	frame1_1 = Frame(frame1, borderwidth=1, relief=RAISED)
	frame1_1.pack(side=LEFT)
	
	w1 = Label(frame1_1, text="Molecules")
	w1.pack(side=TOP)
	
	frame1_1a = Frame(frame1_1)
	frame1_1a.pack(side=TOP)
	
	#List of PyMOL loaded PDB files
	if allNames[0] != "nothing":
		for molecule in allNames:
			radio_button1 = Radiobutton(frame1_1a, text=molecule, value=molecule, variable=v1_name, command=lambda: set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var, p, 1))
			radio_button1.pack(side=TOP, anchor=W)
	#If no loaded PDB files, than add button to choose one
	else:
		w1_1 = Label(frame1_1a, text="Choose PDB file")
		w1_1.pack(side=TOP)
		frame1_1_1 = Frame(frame1_1a)
		frame1_1_1.pack(side=TOP)
		label1 = Label(frame1_1_1, textvariable=v1_name)
		label1.pack(side=LEFT)
		button_e1 = Button(frame1_1_1, text = "Browse", command=lambda: select_file(v1_name))
		button_e1.pack(side=LEFT)
	
	#List of previous projects
	if os.path.isdir(dir_path_dynamics) == True:
		projects = os.listdir(dir_path_dynamics)
	else:
		projects = []
	if projects != ['test_gromacs.txt'] and projects != []:
		w1_2 = Label(frame1_1, text="Previous Projects")
		w1_2.pack(side=TOP)
	
		for molecule in projects:
			if os.path.isdir(dir_path_dynamics + molecule) == True and molecule != "nothing":
				projects1 = []
				molecule1 = molecule.split("_")
				if molecule1[-1] == "multimodel":
					pass
				else:
					molecule1 = molecule.split("-")
					if molecule1[0] == "gromacs":
						pass
					else:
						radio_button1 = Radiobutton(frame1_1, text=molecule, value=molecule, variable=v1_name, command=lambda: set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var, p, 0, check1_button, check_var2))
						radio_button1.pack(side=TOP, anchor=W)
	
	#List of group for final model
	w2 = Label(frame1_1, text="Group")
	w2.pack(side=TOP)
	for group in gromacs.group_list:
		radio_button2 = Radiobutton(frame1_1, text=group[1], value=group[0], variable=v2_group, command=lambda: gromacs2.update(v2_group.get(), gromacs2.box_type, gromacs2.box_distance, gromacs2.box_density))
		radio_button2.pack(side=TOP, anchor=W)
	
	frame1_2 = Frame(frame1, borderwidth=1, relief=RAISED)
	frame1_2.pack(side=LEFT)
	
	#List of available force forces
	w3 = Label(frame1_2, text="Force force", anchor=E)
	w3.pack(side=TOP)

	for force in gromacs.force_list:
		radio_button3 = Radiobutton(frame1_2, text=force[1], value=force[0], variable=v3_force, command=lambda : waterSet(v4_water, water_v, v3_force.get()))
		radio_button3.pack(side=TOP, anchor=W)

	#Label of choosen water model
	w4 = Label(frame1_2, text="Water Model", anchor=E)
	w4.pack(side=TOP)
	
	frame1_2_1 = Frame(frame1_2)
	frame1_2_1.pack(side=TOP)

	#Buttons to choose water model and configure water box
	water_label = Label(frame1_2_1, textvariable=water_v)
	water_label.pack(side=LEFT)
	water_button = Button(frame1_2_1, text = "Choose...", command=lambda : waterChoose(v4_water, water_v, root))
	water_button.pack(side=LEFT)
	waterbox_button = Button(frame1_2_1, text = "Configure", command=lambda: waterBox(root))
	waterbox_button.pack(side=LEFT)
	
	frame1_3 = Frame(frame1)
	frame1_3.pack(side=LEFT)
	
	frame1_3_1 = Frame(frame1_3, borderwidth=1, relief=RAISED)
	frame1_3_1.pack(side=TOP)
	
	w4 = Label(frame1_3_1, text="Configuration")
	w4.pack(side=TOP)
	
	#Button for configuration of MDP files
	em_label = Label(frame1_3_1, text="Energy Minimisation")
	em_label.pack(side=TOP)
	em_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("em", root))
	em_button2.pack(side=TOP)
	if os.path.isfile(dir_path_dynamics + "em.mdp"):
		em_button2.configure(state=DISABLED)
	
	pr_label = Label(frame1_3_1, text="Position Restrained MD")
	pr_label.pack(side=TOP)
	pr_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("pr", root))
	pr_button2.pack(side=TOP)
	if os.path.isfile(dir_path_dynamics + "pr.mdp"):
		pr_button2.configure(state=DISABLED)
	
	md_label = Label(frame1_3_1, text="Molecular Dynamics Simulation")
	md_label.pack(side=TOP)
	md_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("md", root))
	md_button2.pack(side=TOP)
	if os.path.isfile(dir_path_dynamics + "md.mdp"):
		md_button2.configure(state=DISABLED)
	
	#Checkbox and button for position restraints configuration
	check1 = Checkbutton(frame1_3_1, text="Position Restraints", variable=check_var2, command=lambda: RestraintsWindow.check(check_var2.get(), check1_button))
	check1.pack(side=TOP)
	
	check1_button = Button(frame1_3_1, text = "Configure", command=lambda: RestraintsWindow.window(root))
	check1_button.pack(side=TOP)
	check1_button.configure(state=DISABLED)

	frame1_3_2 = Frame(frame1_3, borderwidth=1, relief=RAISED)
	frame1_3_2.pack(side=TOP)
		
	#Jobs configuration button, progress bar and check box
	job_button = Button(frame1_3_2, text = "Jobs to do", command=lambda: jobs_configure(root))
	job_button.pack(side=TOP)
	
	w5 = Label(frame1_3_2, text="Progress Bar")
	w5.pack(side=TOP)
	p = Meter(frame1_3_2, value=0.0)
	p.pack(side=TOP)

	c = Checkbutton(frame1_3_2, text="Molecular Dynamics Simulation from the beginning", variable=check_var, command=lambda: main_status_bar(check_var.get(), p))
	c.pack(side=TOP)
	
	main_status_bar(check_var.get(), p)
	
	frame2 = Frame(root)
	frame2.pack(side=TOP)
	
	#Additional Buttons
	quit_button = Button(frame2, text = "Cancel", command=root.destroy)
	quit_button.pack(side=LEFT)
	
	clean_button = Button(frame2, text = "Clean", command=cleanMessage)
	clean_button.pack(side=LEFT)
	
	help_button = Button(frame2, text = "Help", command=lambda: helpWindow(root))
	help_button.pack(side=LEFT)
	
	save_button = Button(frame2, text = "Save", command=select_file_save)
	save_button.pack(side=LEFT)
	
	load_button = Button(frame2, text = "Load", command=lambda: select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v, check_var))
	load_button.pack(side=LEFT)
	
	count_button = Button(frame2, text = "OK", command=lambda: calculationW.window(root))
	count_button.pack(side=LEFT)
		
	root.mainloop()

##Molecular Dynamics Performing window
class CalculationWindow:
	
	tasks_to_do = 0
	bar_var = ""
	bar_widget = ""
	
	##This function will create main Calculation Window
	def window(self, master):
		print project_name
		master.destroy()
	
		gromacs2.status = ["ok", "Waiting to start"]
	
		root = Tk()
		root.wm_title("Calculation Window")
		frame1 = Frame(root)
		frame1.pack(side=TOP)
		frame2 = Frame(root)
		frame2.pack(side=TOP)
	
		self.bar_var = StringVar(root)
		self.bar_var.set("Ready to start")
	
		w5 = Label(frame1, textvariable=bar_v)
		w5.pack(side=TOP)
		self.bar_widget = Meter(frame1, value=0.0)
		self.bar_widget.pack(side=TOP)
	
		exit_button = Button(frame2, text = "EXIT", command=root.destroy)
		exit_button.pack(side=LEFT)
	
		save_button = Button(frame2, text = "SAVE", command=lambda: select_file_save(1))
		save_button.pack(side=LEFT)
	
		stop_button = Button(frame2, text = "STOP", command=lambda : stop_counting(1))
		stop_button.pack(side=LEFT)
	
		start_button = Button(frame2, text = "START", command=lambda: thread.start_new_thread(dynamics, ()))
		start_button.pack(side=LEFT)
		#Updateing status bar
		tasks_nr = 0.0
		for task in progress.to_do:
			tasks_nr = tasks_nr + task
		self.tasks_to_do = tasks_nr
		thread.start_new_thread(self.bar_update, ())
	
		root.mainloop()

	##This function will update status bar during molecular dynamics simulation
	def bar_update(self):
		global error
		while self.bar_var.get() != "Finished!" and error == "":
			percent = 0.0
			for job in progress.to_do:
				percent = percent + job
			if tasks != 0:
				percent = (tasks - percent) / tasks
			else:
				percent = 1
			self.bar_widget.configure(value=percent)
			self.bar_widget.update_idletasks()
			self.bar_var.set(gromacs2.status[1])
			if stop == 1:
				self.bar_var.set("User Stoped")
			time.sleep(1)
		if self.bar_var.get() == "Finished!":
			print "Finished!"
			self.bar_widget.configure(value=1)
			self.bar_widget.update_idletasks()
			if plugin == 0:
				root = Tk()
				file = tkFileDialog.asksaveasfile(parent=root, mode='w' ,title='Choose final multimodel file to save')
				try:
					os.remove(file.name)
					shutil.copy(dir_path_project + '/' +project_name + "_multimodel.pdb", file.name+".pdb")
				except:
					pass
				root.destroy()
		elif error != "":
			self.bar_var.set("Fatal Error")
			root = Tk()
			root.wm_title("GROMACS Error Message")
			frame = Frame(root)
			frame.pack()
			w = Label(frame, text=error)
			w.pack()
			ok_button = Button(frame, text = "OK", command=root.destroy)
			ok_button.pack()
			root.mainloop()
			error = ""

	##This function will change global value if stop is clicked during simulation
	def stop_counting(self, value):
		global stop
		stop = value

##This function will show current progress in main window if "Jobs from the begining" is unchecked
def main_status_bar(var, bar):
	percent = 0.0
	for job in progress.status:
		percent = percent + job
	percent = percent / 7
	if var == 0:
		bar.configure(value=percent)
		bar.update_idletasks()
		if progress.to_do == [1,1,1,1,1,1,1]:
			progress.to_do_status()
		progress.from_begining = 0
	elif var == 1:
		bar.configure(value=(0.0))
		bar.update_idletasks()
		progress.to_do = [1,1,1,1,1,1,1]
		progress.from_begining = 1

##This function will allow you to choose PDB file if no file is loaded to PyMOL
def select_file():
	root = Tk()
	file = tkFileDialog.askopenfile(parent=root, mode='rb',title='Choose PDB file')
	try:
		name = file.name.split("/")
		name2 = name[-1].split(".")
		##Checking directories
		global project_name, dir_path_project
		project_name = name2[0]
		dir_path_project = dir_path_dynamics + project_name + "/"
		if os.path.isdir(dir_path_project) == False:
			os.makedirs(dir_path_project)
			shutil.copyfile(file.name, dir_path_project + project_name + ".pdb")
			print "pdb_copied"
		set_config_files()	
	except:
		pass
	root.destroy()

def select_file_save(rest_of_work=0):
	if rest_of_work == 1:
		progress.to_do_status()
	root = Tk()
	file = tkFileDialog.asksaveasfile(parent=root, mode='w' ,title='Choose save file')
	if file != None:
		save_file(file.name)	
	root.destroy()

def select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v, check_var):
	root = Tk()
	file = tkFileDialog.askopenfile(parent=root, mode='rb', defaultextension=".tar.bz2" ,title='Choose file to load')
	if file != None:
		new_name = load_file(file.name)
		v1_name.set(new_name)
		v2_group.set(gromacs.group_list[gromacs2.group][0])
		v3_force.set(gromacs.force_list[gromacs2.force-1][0])
		v4_water.set(gromacs.water_list[gromacs2.water-1][0])
		water_v.set(gromacs.water_list[v4_water.get()-1][1])
		check_var.set(progress.from_begining)
		radio_button1 = Radiobutton(frame1_1a, text=new_name, value=new_name, variable=v1_name, command=lambda: set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var, p))
		radio_button1.pack(side=TOP, anchor=W)
	root.destroy()

##Water chooser window
def waterChoose(v4_water, water_v, master):
	root = Toplevel(master)
	root.wm_title("Water Model")
	for water in gromacs.water_list:
		radio_button = Radiobutton(root, text=water[1], value=water[0], variable=v4_water, command = lambda : waterSet(v4_water, water_v))
		radio_button.pack(side=TOP, anchor=W)
	ok_button = Button(root, text = "OK", command=root.destroy)
	ok_button.pack(side=TOP)

##This function will change force field and water model when choosing Force Field in Main Window and also change water model after choosing one in "waterChoose"
def waterSet(v4_water, water_v, force = ""):
	if force == "":
		force = gromacs2.force
	else:
		gromac2.force = force
	gromacs.water_update(force)
	water_v.set(gromacs.water_list[v4_water.get()-1][1])
	gromacs2.water = v4_water.get()
	save_options()

##Water box configuration window
def waterBox(master):
	root1 = Toplevel(master)
	root1.wm_title("Water Box Options")
	root1.wm_geometry("300x250")
	v = StringVar(root1)
	v.set(gromacs2.box_type)
	w = Label(root1, text="Box type")
	w.pack()
	radio_button = Radiobutton(root1, text="triclinic", value="triclinic", variable=v, command = lambda : gromacs2.update(gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="cubic", value="cubic", variable=v, command = lambda : gromacs2.update(gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="dodecahedron", value="dodecahedron", variable=v, command = lambda : gromacs2.update(gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="octahedron", value="octahedron", variable=v, command = lambda : gromacs2.update(gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	w1 = Label(root1, text="Distance")
	w1.pack()
	distance = Entry(root1)
	distance.pack(side=TOP)
	distance.insert(0, gromacs2.box_distance)
	w2 = Label(root1, text="Density [g/L]")
	w2.pack()
	density = Entry(root1)
	density.pack(side=TOP)
	density.insert(0, gromacs2.box_density)
	ok_button = Button(root1, text = "OK", command=lambda : gromacs2.update(gromacs2.group, gromacs2.box_type, distance.get(), density.get(), root1))
	ok_button.pack(side=TOP)

##This is very important function, which sets all necessary configuration files
def set_config_files(name="", v2_group="", v3_force="", v4_water="", water_v="", check_var="", main_bar=0, from_pymol=0, config_button_restraints = "", checkbox_restraints=""):
	global em_file, pr_file, md_file, progress
	if name != "":
		global project_name, dir_path_project
		project_name = name
		dir_path_project = dir_path_dynamics + project_name + '/'
	if os.path.isfile(dir_path_project+"options.pickle") == True:
		load_options()
		v2_group.set(gromacs.group_list[gromacs2.group][0])
		v3_force.set(gromacs.force_list[gromacs2.force-1][0])
		v4_water.set(gromacs.water_list[gromacs2.water-1][0])
		water_v.set(gromacs.water_list[v4_water.get()-1][1])
		if restraintsW.on_off == 0 and config_button_restraints != "":
			checkbox_restraints.set(restraintsW.on_off)
			config_button_restraints.configure(state=DISABLED)
		elif restraintsW.on_off == 1 and config_button_restraints != "":
			checkbox_restraints.set(restraintsW.on_off)
			config_button_restraints.configure(state=ACTIVE)
	else:
		progress = Progress_status(dir_path_project)
		if os.path.isfile(dir_path_dynamics + "em.mdp"):
			shutil.copy(dir_path_dynamics + "em.mdp", dir_path_project + "em.mdp")
			em_file_config = open(dir_path_dynamics + "em.mdp", "r").read()
			em_file = Mdp_config("em.mdp",name, em_file_config, 1)
			print "Found em.mdp file. Using it instead of local configuration."
		else:
			em_file = Mdp_config("em.mdp",em_init_config, 0)
		if os.path.isfile(dir_path_dynamics + "pr.mdp"):
			shutil.copy(dir_path_dynamics + "pr.mdp", dir_path_project + "pr.mdp")
			pr_file_config = open(dir_path_dynamics + "pr.mdp", "r").read()
			pr_file = Mdp_config("pr.mdp",name, pr_file_config, 1)
			print "Found pr.mdp file. Using it instead of local configuration."
		else:
			pr_file = Mdp_config("pr.mdp",pr_init_config, 0)
		if os.path.isfile(dir_path_dynamics + "md.mdp"):
			shutil.copy(dir_path_dynamics + "md.mdp", dir_path_project + "md.mdp")
			md_file_config = open(dir_path_dynamics + "md.mdp", "r").read()
			md_file = Mdp_config("md.mdp",name, md_file_config, 1)
			print "Found md.mdp file. Using it instead of local configuration."
		else:
			md_file = Mdp_config("md.mdp",md_init_config, 0)
	try:
		check_var.set(progress.from_begining)
		if main_bar != 0:
			main_status_bar(check_var.get(), main_bar)
	except:
		pass
	save_options()
	if from_pymol == 1:
		print "cmd saved"
		cmd.save(dir_path_project+project_name+".pdb", project_name) #PyMOL API

#This function will create the window with configuration files based on MDP class
def mdp_configure(config_name, master):
	
	root2 = Toplevel(master)
	
	if config_name == "em":
		options = em_file.options
		root2.wm_title("Energy Minimalization Options")
	elif config_name == "pr":
		options = pr_file.options
		root2.wm_title("Position Restrained MD Options")
	elif config_name == "md":
		options = md_file.options
		root2.wm_title("Molecular Dynamics Simulation Options")

	values_list = []
	for option, value in options:
		frame1 = Frame(root2)
		frame1.pack(side=TOP)
		if option == "emtol":
			l1 = Label(frame1, text="Energy minimizing stuff")
			l1.pack(side=TOP)
		elif option == "Tcoupl":
			l1 = Label(frame1, text="Berendsen temperature and coupling")
			l1.pack(side=TOP)
		elif option == "Pcoupl":
			l1 = Label(frame1, text="Pressure coupling")
			l1.pack(side=TOP)
		elif option == "gen_vel":
			l1 = Label(frame1, text="Generate velocites temperature")
			l1.pack(side=TOP)
		values_list.append(StringVar(root2)) #brilliant
		values_list[-1].set(value)
		if option[0] != ";":
			l = Label(frame1, text=option, width=25, anchor=W)
			l.pack(side=LEFT)
			e = Entry(frame1, textvariable=values_list[-1])
			e.pack(side=LEFT)
	if config_name == "em":
		b = Button(root2, text="OK", command=lambda: mdp_update(values_list, "em", root2))
		b.pack(side=TOP)
	elif config_name == "pr":
		b = Button(root2, text="OK", command=lambda: mdp_update(values_list, "pr", root2))
		b.pack(side=TOP)
	elif config_name == "md":
		b = Button(root2, text="OK", command=lambda: mdp_update(values_list, "md", root2))
		b.pack(side=TOP)

##This function will update MDP class objects alfter closing "mdp_configure" window
def mdp_update(values, mdp, root=""):
	try:
		root.destroy()
	except:
		pass
	values2 = []
	for value in values:
		values2.append(value.get())
	if mdp == "em":
		em_file.update(values2)
	elif mdp == "pr":
		pr_file.update(values2)
	elif mdp == "pr":
		md_file.update(values2)

##This function will create Jobs configuration window
def jobs_configure(master):
	root5 = Toplevel(master)
	root5.wm_title("Jobs configuration")
	check_var1 = IntVar(root5)
	check_var1.set(progress.to_do[0])
	check_var2 = IntVar(root5)
	check_var2.set(progress.to_do[1])
	v1 = IntVar(root5)
	v1.set(1)
	check_var3 = IntVar(root5)
	check_var3.set(progress.to_do[2])
	check_var4 = IntVar(root5)
	check_var4.set(progress.to_do[3])
	check_var5 = IntVar(root5)
	check_var5.set(progress.to_do[4])
	check_var6 = IntVar(root5)
	check_var6.set(progress.to_do[5])
	check_var7 = IntVar(root5)
	check_var7.set(progress.to_do[6])
	
	frame1 = Frame(root5)
	frame1.pack(side=TOP)

	c1 = Checkbutton(frame1, text="Save configuration files", variable=check_var1, command=lambda: progress.to_do_update(0, check_var1.get()))
	c1.pack(side=TOP, anchor=W)
	
	c2 = Checkbutton(frame1, text="Generate topology file from pdb", variable=check_var2, command=lambda: progress.to_do_update(1, check_var2.get()))
	c2.pack(side=TOP, anchor=W)
	
	r1 = Radiobutton(frame1, text="Use pdb2gmx tool", value=1, variable=v1, command=lambda: test(v1.get()))
	r1.pack(side=TOP, anchor=W)
	#Not yet available
	r2 = Radiobutton(frame1, text="Use x2top tool [coming soon]", value=2, variable=v1, command=lambda: test(v1.get()))
	r2.pack(side=TOP, anchor=W)
	r2.configure(state=DISABLED)
	
	c3 = Checkbutton(frame1, text="Adding Water Box", variable=check_var3, command=lambda: progress.to_do_update(2, check_var3.get()))
	c3.pack(side=TOP, anchor=W)
	
	c4 = Checkbutton(frame1, text="Energy minimalization", variable=check_var4, command=lambda: progress.to_do_update(3, check_var4.get()))
	c4.pack(side=TOP, anchor=W)
	
	c5 = Checkbutton(frame1, text="Position Restrained MD", variable=check_var5, command=lambda: progress.to_do_update(4, check_var5.get()))
	c5.pack(side=TOP, anchor=W)
	
	c6 = Checkbutton(frame1, text="Molecular Dynamics Simulation", variable=check_var6, command=lambda: progress.to_do_update(5, check_var6.get()))
	c6.pack(side=TOP, anchor=W)
	
	c7 = Checkbutton(frame1, text="Generate multimodel PDB", variable=check_var7, command=lambda: progress.to_do_update(6, check_var7.get()))
	c7.pack(side=TOP, anchor=W)
	
	b1 = Button(root5, text="OK", command=root5.destroy)
	b1.pack(side=TOP)
	
##Restraints
##This class is resposible for graphic edition of restraints
class RestraintsWindow:
	
	atom_list = []
	check_var = ""
	on_off = 0
	
	def __init__(self):
		print "Need to prepare here some atome list"
	
	##This function will create main window for restraints
	def window(self, master):
		root1 = Toplevel(master)
		root1.wm_title("Restraints Configure")
	
		sw = ScrolledWindow(root1, scrollbar=Y) # just the vertical scrollbar
		sw.pack(fill=BOTH, expand=1)
	
		self.check_var = IntVar(sw.window)
		self.check_var.set(gromacs2.restraints_nr)
	
		number = 0
		for group in gromacs.restraints:
			select = Radiobutton(sw.window, text=group[0], value=number, variable=self.check_var)
			select.pack()
			text = Text(sw.window)
			text.insert(END, group[1])
			text.pack()
			self.atom_list.append(text)
			number = number + 1
	
		if plugin == 1:
			select1 = Radiobutton(sw.window, text="[ PyMol Selected ]", value=number, variable=self.check_var)
			select1.pack()
			text1 = Text(sw.window)
		
			stored.list=[]
			cmd.iterate("(sele)","stored.list.append(ID)") #PyMOL API
		
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
		
		ok_button = Button(sw.window, text="OK", command=lambda : self.index(root1))
		ok_button.pack()

	##This function will modyfie index_dynamics.ndx file based on user choosed restraints		
	def index(self, root_to_kill=""):
		index_nr = self.check_var.get()
		gromacs2.restraints_nr = index_nr
		text = self.atom_list[index_nr]
		if index_nr < len(gromacs.restraints):
			gromacs.restraints[index_nr][1] = text.get(1.0, END)
		index_file = open("index_dynamics.ndx", "w")
		index_file.write("[ Dynamics Selected ]\n"+text.get(1.0, END))
		index_file.close()
		if root_to_kill != "":
			root_to_kill.destroy()
			
	##This function will activ or disable restraints button in main window based on check box
	def check(self, check, config_button):
		global restraintsW.on_off
		if check == 1:
			config_button.configure(state=ACTIVE)
			md_file.options[2][0] = "define"
			md_file.update()
			self.on_off = 1
			gromacs.restraintsW.on_offindex()
		elif check == 0:
			config_button.configure(state=DISABLED)
			md_file.options[2][0] = ";define"
			md_file.update()
			self.on_off = 0

##Help window
def helpWindow(master):
	root = Toplevel(master)
	root.wm_title("Help Window")
	frame = Frame(root)
	frame.pack()
		
	w = Label(frame, text=help_option())
	w.pack()
	ok_button = Button(frame, text = "OK", command=root.destroy)
	ok_button.pack()

##Clean message in tkMessageBox
def cleanMessage():
	tkMessageBox.showinfo("Clean", "Temporary files are now removed!\nPlease restart plugin.")
	clean_option()

##--Comand Line Interface--
def dynamics(help_clean = ""):
	file_path = dir_path_project + project_name
	force = str(gromacs2.force) + "\n" + str(gromacs2.water)
	group = str(gromacs2.group)
	global status, stop
	stop = 0
	print "Starting PyMOL plugin 'dynamics' ver."+plugin_ver+" by Tomasz Makarewicz"
	
	##If help is called
	if help_name.count(help_clean) == 1:
		print help_option()
		status = ["fail", "Help printed"]
	##If clean is called
	elif clean_name.count(help_clean) == 1:
		clean_option()
		status = ["fail", "Cleaned"]
	
	##Starting real calculations
	##Checking GROMACS
	if status[0] == "ok":
		status = gromacs.status
	
	##Saving configuration files
	if status[0] == "ok" and stop == 0 and progress.to_do[0] == 1:
		mdp_files()
		if status[0] == "ok":
			progress.status[0] = 1
			progress.to_do[0] = 0
			save_options()

	##Checking variables - temporary disabled
	#if status[0] == "ok":
	#	status = check_variable(name, force1, force2, group)
	#elif status[0] == "fail":
	#	pass
	
	##Counting topology
	if status[0] == "ok" and stop == 0 and progress.to_do[1] == 1:
		gromacs2.pdb2top(file_path, force)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[1] = 1
			progress.to_do[1] = 0
			save_options()

	##Adding water box
	if status[0] == "ok" and stop == 0 and progress.to_do[2] == 1:
		gromacs2.waterbox(file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[2] = 1
			progress.to_do[2] = 0
			save_options()
	
	##EM	
	if status[0] == "ok" and stop == 0 and progress.to_do[3] == 1:
		gromacs2.em(file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[3] = 1
			progress.to_do[3] = 0
			save_options()
	
	##PR
	if status[0] == "ok" and stop == 0 and progress.to_do[4] == 1:
		gromacs2.pr(file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[4] = 1
			progress.to_do[4] = 0
			save_options()
	
	##Restraints
	if status[0] == "ok" and stop == 0 and restraintsW.on_off == 1 and progress.to_do[5] == 1:
		gromacs2.restraints()
		status = gromacs2.status
	
	##MD
	if status[0] == "ok" and stop == 0 and progress.to_do[5] == 1:
		gromacs2.md(file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[5] = 1
			progress.to_do[5] = 0
			save_options()
	
	##Trjconv
	if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
		gromacs2.trjconv(file_path, group)
		status = gromacs2.status
	
	##Showing multimodel
	if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
		show_multipdb()
		progress.status[6] = 1
		progress.to_do[6] = 0
		save_options()
	elif status[0] == "fail":
		print status[1]
		if help_name.count(help_clean) != 1 and clean_name.count(help_clean) != 1:
			error_message()
		status = ["ok", ""]
	#This is for runnig program outside PyMOL
	if plugin == 0:
		return dir_path_project

##Checking if given varaibles are correct - depreciated
def check_variable(force1, force2, group):
	status = ["ok", ""]
	try:
		int(force1)
	except:
		status = ["fail", "Wrong Force force"]
	try:
		int(force2)
	except:
		status = ["fail", "Wrong Water Model"]
	try:
		int(group)
	except:
		status = ["fail", "Wrong Group"]
	return status

##Saving configuration files
def mdp_files():
	em_file.save_file()
	pr_file.save_file()
	md_file.save_file()

##Show multimodel PDB file in PyMOOL
def show_multipdb():
	if plugin == 1:
		try:
			cmd.hide("everything", project_name) #PyMOL API
		except:
			pass
		cmd.load(project_name+"_multimodel.pdb") #PyMOL API

##Saving tar.bz file
def save_file(destination_path):
	print "Saving"
	import tarfile
	save_options()
	tar = tarfile.open(destination_path+".tar.bz2", "w:bz2")
	tar.add(dir_path_project, recursive=True, arcname=name[-2])
	tar.close()
	os.remove(destination_path)

##Load tar.bz file
def load_file(file_path):
	print "Loading"
	import tarfile
	tar = tarfile.open(file_path, "r:bz2")
	names = tar.getnames()
	#Backup same name folder if file is loaded
	if os.path.isdir(dir_path_dynamics + names[0]) == True:
		back_folder = dir_path_dynamics + names[0] + "_back"
		while os.path.isdir(back_folder) == True:
			back_folder = back_folder + "_b"
		os.rename(dir_path_dynamics + names[0], back_folder)
	tar.extractall(dir_path_dynamics)
	global dir_path_project, project_name
	project_name = names[0]
	dir_path_project = dir_path_dynamics + project_name + "/"
	load_options()

##Save all settings to options.pickle file
def save_options():
	if os.path.isdir(dir_path_project) == False:
		os.makedirs(dir_path_project)
	destination_option = file(dir_path_project + "options.pickle", "w")
	pickle.dump([plugin_ver, gromacs.version, gromacs2, em_file, pr_file, md_file, progress, restraintsW], destination_option)
	del destination_option

##Load all settings from options.pickle file	
def load_options():
	global gromacs, gromacs2, em_file, pr_file, md_file, progress, restraintsW
	
	pickle_file = file(dir_path_project +"options.pickle")
	options = pickle.load(pickle_file)
	
	print "Project was created for Dynamics PyMOL Plugin"+options[0]+" and "+options[1]
	if gromacs.version != options[1]:
		print "GROMACS versions not identical."
	if options[0][1:4] == "1.0":
		print "1.0 compatibility layer"
		gromacs2 = options[2]
		em_file = options[3]
		pr_file = options[4]
		md_file = options[5]
		progress = options[6]
		restraintsW.on_off = options[7]

	elif options[0][1:4] == "1.1":
		print "1.1 compatibility layer"
		gromacs2 = options[2]
		em_file = options[3]
		pr_file = options[4]
		md_file = options[5]
		progress = options[6]
		restraintsW = options[7]

##Text for "Help"
def help_option():
	help_message = """This is the dynamics PyMOL Plugin.
This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
Contributors:
- Tomasz Makarewicz (tomaszm@biotech.ug.edu.pl)

Full manual is available to you on project website: ...
or as a file: /usr/share/doc/dynamics-pymol-plugin/manual.pdf

The purpose of this plugin is to perform molecular dynamics simulation by GROMACS using easy graphical tool and powerful molecular viewer.

To use this program run it as a PyMOL plugin.
Choose molecule (PDB) for which you want to perform molecular dynamics simulation (left column).
Choose force force and water model options in the middle column.
Choose any additional options in the right column.
Press OK button.
Click Start button and wait till calculations are finished.
You will see animated model in PyMOL viewer."""
	return help_message

##Clean function
def clean_option():
	shutil.rmtree(dir_path_dynamics)
	print "Temporary files are now removed."

##If molecular dynamics simulation fails, this function will show the error
def error_message():
	
	global error
	
	log = open(dir_path_project+"log.txt","r")
	log_list = log.readlines()
	
	error_start_line = 0
	while log_list[error_start_line] != "Fatal error:\n":
		error_start_line = error_start_line + 1
	error_end_line = error_start_line
	while log_list[error_end_line] != "-------------------------------------------------------\n":
		error_end_line = error_end_line + 1
	error_list = log_list[error_start_line:error_end_line]
	error = ""
	for line in error_list:
		error = error + line
	print error
	
##--PyMOL Shell Interface-- - depreciated
if plugin == 1:
	cmd.extend("dynamics", dynamics) #PyMOL API
