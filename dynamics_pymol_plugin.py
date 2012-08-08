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
plugin_ver = " 1.0.3pre"

stop = 0
restraints_var = 0
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

dir_path_dynamics = os.getenv("HOME")+'/.dynamics/'
dir_path_project = dir_path_dynamics+'nothing'
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
	force_field_list = []
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

		print "Reading available force fields"	
		force_field_start_line = 0
		while lista_gromacs[force_field_start_line] != "Select the Force Field:\n":
			force_field_start_line = force_field_start_line + 1
		force_field_start_line = force_field_start_line + 2
		force_field_end_line = force_field_start_line
		while lista_gromacs[force_field_end_line] != "\n":
			force_field_end_line = force_field_end_line + 1
		force_field_list = lista_gromacs[force_field_start_line:force_field_end_line]
		force_field_list2 = []
		number = 1
		for force in force_field_list:
			force_field_list2.append([number, force[:-1]])
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
		
		self.force_field_list = force_field_list2
		self.water_list = water_list2
		self.group_list = group_list2
	
	##This function will update water list if force field is changed.
	def water_update(self, force_field_number):
		print "Updateing available water models"
		subprocess.call(self.path+"echo -e '"+str(force_field_number)+"\n1' | pdb2gmx &> "+dir_path_dynamics+"test_gromacs.txt", executable="/bin/bash", shell=True)
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
		if os.path.isdir(dir_path_project) == True:
			save1()
		return water_list2
	
	##This function will read atoms group for restraints for current molecule.	
	def restraints_index(self):
		self.restraints = []
		if os.path.isdir(dir_path_project) == False:
			os.makedirs(dir_path_project)
		os.chdir(dir_path_project)
		name = dir_path_project.split("/")
		name = name[-2]
		subprocess.call(self.path+"echo q | make_ndx -f "+name+".pdb -o index.ndx &> restraints.log", executable="/bin/bash", shell=True)	
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
	field = 1
	water = 1
	group = 1
	box_type = "triclinic"
	box_distance = "0.5"
	box_density = "1000"
	restraints_nr = 1

	##This function will change given variabless stored by the class (unnecesary - only for backward compatybility, need to be removed)
	def update(self, field, water, group, box_type, box_distance, box_density, root=""):
		#Close mother window if present
		try:
			root.destroy()
		except:
			pass
		
		self.field = field
		self.water = water
		self.group = group
		self.box_type = box_type
		self.box_distance = box_distance
		self.box_density = box_density
		save1()
		print "gromacs update"
	
	##This function will create initial topology and triectory using pdb file and choosen force field
	def pdb2top(self, name, file_path, field):
		status = ["ok", "Calculating topology using Force Fields"]
		self.status = status
		try:
			os.remove(name+".gro")
			os.remove(name+".top")
		except:
			pass
		print "Calculating topology using Force Fields"
		Pdb2gmx = subprocess.call(gromacs.path+"echo -e '"+field+"' | pdb2gmx -f "+name+".pdb -o "+name+".gro -p "+name+".top &> log.txt",
		executable="/bin/bash", shell=True)

		if os.path.isfile(file_path+".gro") == True:
			status = ["ok", ""]
		else:
			status = ["fail", "Warning. Trying to ignore unnecessary hydrogen atoms."]
			print status[1]
			Pdb2gmx = subprocess.call(gromacs.path+"echo -e '"+field+"' | pdb2gmx -ignh -f "+name+".pdb -o "+name+".gro -p "+name+".top &> log.txt",
			executable="/bin/bash", shell=True)

		if os.path.isfile(file_path+".gro") == True:
			status = ["ok", "Calculated topology using Force Fields"]
		else:
			status = ["fail", "Force field unable to create topology file"]
		
		self.status = status
	
	##This function will create and add waterbox.
	def waterbox(self, name, file_path):
		status = ["ok", "Adding Water Box"]
		self.status = status
		box_type = "-bt "+self.box_type+" "
		distance = "-d "+self.box_distance+" "
		density = "-density "+self.box_density
		
		try:
			os.remove(name+"1.gro")
		except:
			pass
		
		print "Generating water_box"
		Editconf = subprocess.call(gromacs.path+"editconf -f "+name+".gro -o "+name+"1.gro "+box_type+distance+density+" &>> log.txt",
		executable="/bin/bash", shell=True)

		print "Adding Water Box"
		Genbox = subprocess.call(gromacs.path+"genbox -cp "+name+"1.gro -cs -o "+name+"_b4em.gro -p "+name+".top &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"1.gro") == True:
			status = ["ok", "Added Water Box"]
		else:
			status = ["fail", "Unable to add waterbox"]
		self.status = status
	
	##This function will perform energy minimalization	
	def em(self, name, file_path):
		status = ["ok", "Energy Minimalization"]
		self.status = status
		
		try:
			os.remove(name+"_em.tpr")
			os.remove(name+"_em.trr")
		except:
			pass

		print "Energy Minimalization"		
		Grompp = subprocess.call(gromacs.path+"grompp -f em -c "+name+"_b4em -p "+name+" -o "+name+"_em &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+name+"_em -o "+name+"_em -c "+name+"_b4pr -v &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"_em.tpr") == True:
			status = ["ok", "Energy Minimalized"]
		else:
			status = ["fail", "Unable to perform Energy Minimalization"]
		self.status = status
	
	##This function will perform position restrained MD
	def pr(self, name, file_path):
		status = ["ok", "Position Restrained MD"]
		self.status = status
		
		try:
			os.remove(name+"_pr.tpr")
			os.remove(name+"_pr.trr")
		except:
			pass
		
		print "Position Restrained MD"
		Grompp = subprocess.call(gromacs.path+"grompp -f pr -c "+name+"_b4pr -r "+name+"_b4pr -p "+name+" -o "+name+"_pr &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+name+"_pr -o "+name+"_pr -c "+name+"_b4md -v &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile(file_path+"_pr.tpr") == True:
			status = ["ok", "Position Restrained MD finished"]
		else:
			status = ["fail", "Unable to perform Position Restrained"]
		self.status = status
	
	##This function will create posre.itp file for molecular dynamics simulation with choosen atoms if restraints were selected
	def restraints(self, name):
		status = ["ok", "Adding Restraints"]
		self.status = status
		
		try:
			os.remove("posre_2.itp")
		except:
			pass
			
		print "Adding Restraints"
		Genrestr = subprocess.call("echo 0 | genrestr -f "+name+".pdb -o posre_2.itp -n index_dynamics.ndx &>> log.txt", executable="/bin/bash", shell=True)
		
		if os.path.isfile("posre_2.itp") == True:
			status = ["ok", "Added Restraints"]
			os.remove("posre.itp")
			shutil.copy("posre_2.itp", "posre.itp")
		else:
			status = ["fail", "Unable to create restraints file"]
		self.status = status
	
	##This function will perform position final molecular dynamics simulation
	def md(self, name, file_path):
		status = ["ok", "Molecular Dynamics Simulation"]
		self.status = status
		
		try:
			os.remove(name+"_md.tpr")
			os.remove(name+"_md.trr")
		except:
			pass
		
		print "Molecular Dynamics Simulation"
		Grompp = subprocess.call(gromacs.path+"grompp -f md -c "+name+"_b4md  -p "+name+" -o "+name+"_md &>> log.txt", executable="/bin/bash", shell=True)

		Mdrun = subprocess.call(gromacs.path+"mdrun -nice 4 -s "+name+"_md -o "+name+"_md -c "+name+"_after_md -v &>> log.txt", executable="/bin/bash", shell=True)
	
		if os.path.isfile(file_path+"_md.tpr") == True:
			status = ["ok", "Molecular Dynamics Simulation finished"]
		else:
			status = ["fail", "Unable to perform Molecular Dynamics Simulation"]
		self.status = status
	
	##This function will convert final results to multimodel pdb file
	def trjconv(self, name, file_path, group):
		status = ["ok", "Creating Multimodel PDB"]
		self.status = status
		
		try:
			os.remove(name+"_multimodel.pdb")
		except:
			pass
		if os.path.isfile(name+"_multimodel.pdb") == True:
			os.remove(name+"_multimodel.pdb")
		
		print "Creating Multimodel PDB"
		Trjconv = subprocess.call(gromacs.path+"echo "+group+" | trjconv -f "+name+"_md.trr -s "+name+"_md.tpr -app -o "+name+"_multimodel.pdb &>> log.txt", executable="/bin/bash", shell=True)
		
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
	
	def __init__(self, file_name, name, init_config, external_file=0):
		if external_file == 0:
			self.config = """title = """+name+"_"+file_name+"\n"+init_config
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
		save1()
	
	def save_file(self, dir_path):
		mdp = open(dir_path+self.file_name, "w")
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

##init function - puts plugin into menu and starts it () after clicking.
def __init__(self):
	self.menuBar.addmenuitem("Plugin", "command", "dynamics"+plugin_ver, label = "dynamics"+plugin_ver,
	command = lambda s=self: dynamicsDialog(s))
	
##Master menu window
class MasterWindow:
	
	def __init__(self, master):
		
		##Detect list of PyMOL loaded PDB files if no files than list "nothing"
		if plugin == 1:
			allNames = cmd.get_names("objects") #PyMOL API
			e1 = StringVar(master)
			e1.set("")
			if allNames == []:
				allNames = ["nothing"]
			else:
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
		moleculeName = allNames[0]
		v1_name = StringVar(master)
		v1_name.set(moleculeName)
		
		groupName = gromacs.group_list[1][0]
		v2_group = IntVar(master)
		v2_group.set(groupName)
		
		force_fieldName = gromacs.force_field_list[0][0]
		v3_force = IntVar(master)
		v3_force.set(force_fieldName)
		
		waterName = gromacs.water_list[0][0]
		v4_water = IntVar(master)
		v4_water.set(waterName)

		water_v = StringVar(master)
		water_v.set(gromacs.water_list[v4_water.get()-1][1])
		
		check_var = IntVar(master)
		
		check_var2 = IntVar(master)
		check_var2.set(restraints_var)
		
		#Initial configuration
		set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var)
		
		##Start drawing interface
		frame0 = Frame(master)
		frame0.pack(side=TOP)
		
		w_version = Label(frame0, text=gromacs.version)
		w_version.pack(side=TOP)
		
		frame1 = Frame(master)
		frame1.pack(side=TOP)
		
		frame1_1 = Frame(frame1, borderwidth=1, relief=RAISED)
		frame1_1.pack(side=LEFT)
		
		w1 = Label(frame1_1, text="Molecules")
		w1.pack(side=TOP)
		
		global dir_path_project
		dir_path_project = dir_path_dynamics + v1_name.get().lower() + "/"
		
		frame1_1a = Frame(frame1_1)
		frame1_1a.pack(side=TOP)
		
		#List of PyMOL loaded PDB files
		if allNames[0] != "nothing":
			for molecule in allNames:
				radio_button1 = Radiobutton(frame1_1a, text=molecule, value=molecule, variable=v1_name, command=lambda: set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var, p, 1))
				radio_button1.pack(side=TOP, anchor=W)
		#If no loaded PDB files, than add button to choose one
		else:
			#Change entry to label
			w1_1 = Label(frame1_1a, text="Choose PDB file")
			w1_1.pack(side=TOP)
			frame1_1_1 = Frame(frame1_1a)
			frame1_1_1.pack(side=TOP)
			e1 = Entry(frame1_1_1)
			e1.pack(side=LEFT)
			button_e1 = Button(frame1_1_1, text = "Browse", command=lambda: select_file(e1, v1_name))
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
			radio_button2 = Radiobutton(frame1_1, text=group[1], value=group[0], variable=v2_group, command=lambda: gromacs2.update(v3_force.get(), v4_water.get(), v2_group.get(), gromacs2.box_type, gromacs2.box_distance, gromacs2.box_density))
			radio_button2.pack(side=TOP, anchor=W)
		
		frame1_2 = Frame(frame1, borderwidth=1, relief=RAISED)
		frame1_2.pack(side=LEFT)
		
		#List of available force fields
		w3 = Label(frame1_2, text="Force Field", anchor=E)
		w3.pack(side=TOP)

		for force_field in gromacs.force_field_list:
			radio_button3 = Radiobutton(frame1_2, text=force_field[1], value=force_field[0], variable=v3_force, command=lambda : waterSet(v4_water, water_v, 1, v3_force.get()))
			radio_button3.pack(side=TOP, anchor=W)

		#Label of choosen water model
		w4 = Label(frame1_2, text="Water Model", anchor=E)
		w4.pack(side=TOP)
		
		frame1_2_1 = Frame(frame1_2)
		frame1_2_1.pack(side=TOP)

		#Buttons to choose water model and configure water box
		water_label = Label(frame1_2_1, textvariable=water_v)
		water_label.pack(side=LEFT)
		water_button = Button(frame1_2_1, text = "Choose...", command=lambda : waterChoose(v4_water, water_v, master)) #!!!
		water_button.pack(side=LEFT)
		waterbox_button = Button(frame1_2_1, text = "Configure", command=lambda: waterBox(master))
		waterbox_button.pack(side=LEFT)
		
		frame1_3 = Frame(frame1)
		frame1_3.pack(side=LEFT)
		
		frame1_3_1 = Frame(frame1_3, borderwidth=1, relief=RAISED)
		frame1_3_1.pack(side=TOP)
		
		w4 = Label(frame1_3_1, text="Configuration")
		w4.pack(side=TOP)
		
		#Buttonf for configuration of MDP files
		em_label = Label(frame1_3_1, text="Energy Minimisation")
		em_label.pack(side=TOP)
		em_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("em", master))
		em_button2.pack(side=TOP)
		if os.path.isfile(dir_path_dynamics + "em.mdp"):
			em_button2.configure(state=DISABLED)
		
		pr_label = Label(frame1_3_1, text="Position Restrained MD")
		pr_label.pack(side=TOP)
		pr_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("pr", master))
		pr_button2.pack(side=TOP)
		if os.path.isfile(dir_path_dynamics + "pr.mdp"):
			pr_button2.configure(state=DISABLED)
		
		md_label = Label(frame1_3_1, text="Molecular Dynamics Simulation")
		md_label.pack(side=TOP)
		md_button2 = Button(frame1_3_1, text = "Configure", command=lambda: mdp_configure("md", master))
		md_button2.pack(side=TOP)
		if os.path.isfile(dir_path_dynamics + "md.mdp"):
			md_button2.configure(state=DISABLED)
		
		#Check-box and button for position restraints configuration
		check1 = Checkbutton(frame1_3_1, text="Position Restraints", variable=check_var2, command=lambda: restraints(check_var2.get(), check1_button))
		check1.pack(side=TOP)
		
		check1_button = Button(frame1_3_1, text = "Configure", command=lambda: restraints_window(master))
		check1_button.pack(side=TOP)
		check1_button.configure(state=DISABLED)

		frame1_3_2 = Frame(frame1_3, borderwidth=1, relief=RAISED)
		frame1_3_2.pack(side=TOP)
		
		#Jobs configuration button, progress bar and check box
		job_button = Button(frame1_3_2, text = "Jobs to do", command=lambda: jobs_configure(master))
		job_button.pack(side=TOP)
		
		w5 = Label(frame1_3_2, text="Progress Bar")
		w5.pack(side=TOP)
		p = Meter(frame1_3_2, value=0.0)
		p.pack(side=TOP)
		
		#Set configure if PDB file was loaded in PyMOL
		if allNames[0] != "nothing":
			save1()
			set_config_files(v1_name.get(), v2_group, v3_force, v4_water, water_v, check_var, p, 1)

		c = Checkbutton(frame1_3_2, text="Molecular Dynamics Simulation from the beginning", variable=check_var, command=lambda: main_status_bar(check_var.get(), p))
		c.pack(side=TOP)
		
		main_status_bar(check_var.get(), p)
		
		frame2 = Frame(master)
		frame2.pack(side=TOP)
		
		#Additional Buttons
		quit_button = Button(frame2, text = "Cancel", command=master.destroy)
		quit_button.pack(side=LEFT)
		
		clean_button = Button(frame2, text = "Clean", command=cleanMessage)
		clean_button.pack(side=LEFT)
		
		help_button = Button(frame2, text = "Help", command=lambda: helpWindow(master))
		help_button.pack(side=LEFT)
		
		save_button = Button(frame2, text = "Save", command=lambda: select_file_save(v1_name.get()))
		save_button.pack(side=LEFT)
		
		load_button = Button(frame2, text = "Load", command=lambda: select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v, check_var))
		load_button.pack(side=LEFT)
	
		count_button = Button(frame2, text = "OK", command=lambda: startMessage(v1_name.get(), v2_group.get(), v3_force.get(), v4_water.get(), master))
		count_button.pack(side=LEFT)

##--Graphic Interface--
def dynamicsDialog(app):
	root = Tk()
	root.wm_title("Dynamics"+plugin_ver)
	app = MasterWindow(root)
	root.mainloop()

##Molecular Dynamics Performing window
def startMessage(moleculeName, group_nr, force_field_nr, water_nr, master):
	print moleculeName
	master.destroy()
	
	gromacs2.group = group_nr
	gromacs2.field = force_field_nr
	gromacs2.water = water_nr
	gromacs2.status = ["ok", "Waiting to start"]
	
	root = Tk()
	root.wm_title("Counting")
	frame1 = Frame(root)
	frame1.pack(side=TOP)
	frame2 = Frame(root)
	frame2.pack(side=TOP)
	
	bar_v = StringVar(root)
	bar_v.set("Ready to start")
	
	w5 = Label(frame1, textvariable=bar_v)
	w5.pack(side=TOP)
	p = Meter(frame1, value=0.0)
	p.pack(side=TOP)
	
	exit_button = Button(frame2, text = "EXIT", command=root.destroy)
	exit_button.pack(side=LEFT)
	
	save_button = Button(frame2, text = "SAVE", command=lambda: select_file_save(moleculeName, 1))
	save_button.pack(side=LEFT)
	
	stop_button = Button(frame2, text = "STOP", command=lambda : stop_counting(1))
	stop_button.pack(side=LEFT)
	
	start_button = Button(frame2, text = "START", command=lambda: thread.start_new_thread(dynamics, (moleculeName,)))
	start_button.pack(side=LEFT)
	#Updateing status bar
	tasks_nr = 0.0
	for task in progress.to_do:
		tasks_nr = tasks_nr + task
	thread.start_new_thread(bar_update, (bar_v, p, tasks_nr, moleculeName))
	
	root.mainloop()

##This function will update status bar during molecular dynamics simulation
def bar_update(bar_var, bar_widget, tasks, name):
	global error
	while bar_var.get() != "Finished!" and error == "":
		percent = 0.0
		for job in progress.to_do:
			percent = percent + job
		if tasks != 0:
			percent = (tasks - percent) / tasks
		else:
			percent = 1
		bar_widget.configure(value=percent)
		bar_widget.update_idletasks()
		bar_var.set(gromacs2.status[1])
		if stop == 1:
			bar_var.set("User Stoped")
		time.sleep(1)
	if bar_var.get() == "Finished!":
		print "Finished!"
		bar_widget.configure(value=1)
		bar_widget.update_idletasks()
		if plugin == 0:
			root = Tk()
			file = tkFileDialog.asksaveasfile(parent=root, mode='w' ,title='Choose final multimodel file to save')
			try:
				os.remove(file.name)
				shutil.copy(dir_path_project + name + "_multimodel.pdb", file.name+".pdb")
			except:
				pass
			root.destroy()
	elif error != "":
		bar_var.set("Fatal Error")
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

##This function will change global value if stop is clicked during simulation
def stop_counting(value):
	global stop
	stop = value

##This function will allow you to choose PDB file if no file is loaded to PyMOL
def select_file(entry, project_name):
	root = Tk()
	file = tkFileDialog.askopenfile(parent=root, mode='rb',title='Choose PDB file')
	try:
		entry.insert(0, file.name)
		entry.delete(0, END)
		entry.insert(0, file.name)
		name = file.name.split("/")
		name2 = name[-1].split(".")
		project_name.set(name2[0])
		##Checking directories
		global dir_path_project
		dir_path_project = dir_path_dynamics + name2[0].lower() + "/"
		if os.path.isdir(dir_path_project) == False:
			os.makedirs(dir_path_project)
			shutil.copyfile(file.name, dir_path_project + name2[0].lower() + ".pdb")
			print "pdb_copied"
		set_config_files(name2[0])	
	except:
		pass
	root.destroy()

def select_file_save(moleculeName, rest_of_work=0):
	if rest_of_work == 1:
		progress.to_do_status()
	root = Tk()
	file = tkFileDialog.asksaveasfile(parent=root, mode='w' ,title='Choose save file')
	if file != None:
		save(file.name)	
	root.destroy()

def select_file_load(frame1_1a, v1_name, v2_group, v3_force, v4_water, water_v, check_var):
	root = Tk()
	file = tkFileDialog.askopenfile(parent=root, mode='rb', defaultextension=".tar.bz2" ,title='Choose file to load')
	if file != None:
		new_name = load(file.name)
		v1_name.set(new_name)
		v2_group.set(gromacs.group_list[gromacs2.group][0])
		v3_force.set(gromacs.force_field_list[gromacs2.field-1][0])
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
	waterNr = gromacs.water_list[0][0]
	v_nr = IntVar(root)
	v_nr.set(waterNr)
	for water in gromacs.water_list:
		radio_button = Radiobutton(root, text=water[1], value=water[0], variable=v_nr, command = lambda : waterSet(v4_water, water_v, v_nr.get()))
		radio_button.pack(side=TOP, anchor=W)
	ok_button = Button(root, text = "OK", command=root.destroy)
	ok_button.pack(side=TOP)

##This function will change water model after choosing one in "waterChoose"
def waterSet(v4_water, water_v, water_nr, force = ""):
	if force == "":
		force = gromacs2.field
	gromacs.water_update(force)
	v4_water.set(water_nr)
	water_v.set(gromacs.water_list[water_nr-1][1])
	gromacs2.update(force, water_nr, gromacs2.group, gromacs2.box_type, gromacs2.box_distance, gromacs2.box_density)
	save1()

##Water box configuration window
def waterBox(master):
	root1 = Toplevel(master)
	root1.wm_title("Water Box Options")
	root1.wm_geometry("300x250")
	v = StringVar(root1)
	v.set(gromacs2.box_type)
	w = Label(root1, text="Box type")
	w.pack()
	radio_button = Radiobutton(root1, text="triclinic", value="triclinic", variable=v, command = lambda : gromacs2.update(gromacs2.field, gromacs2.water, gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="cubic", value="cubic", variable=v, command = lambda : gromacs2.update(gromacs2.field, gromacs2.water, gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="dodecahedron", value="dodecahedron", variable=v, command = lambda : gromacs2.update(gromacs2.field, gromacs2.water, gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
	radio_button.pack(side=TOP, anchor=W)
	radio_button = Radiobutton(root1, text="octahedron", value="octahedron", variable=v, command = lambda : gromacs2.update(gromacs2.field, gromacs2.water, gromacs2.group, v.get(), gromacs2.box_distance, gromacs2.box_density))
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
	ok_button = Button(root1, text = "OK", command=lambda : gromacs2.update(gromacs2.field, gromacs2.water, gromacs2.group, gromacs2.box_type, distance.get(), density.get(), root1))
	ok_button.pack(side=TOP)

##This is very important function, which sets all necessary configuration files
def set_config_files(name, v2_group="", v3_field="", v4_water="", water_v="", check_var="", main_bar=0, from_pymol=0, config_button_restraints = "", checkbox_restraints=""):
	global em_file, pr_file, md_file, progress, dir_path_project
	name = name.lower()
	dir_path_project = dir_path_dynamics + name + "/"

	if os.path.isfile(dir_path_project+"options.pickle") == True:
		load1()
		v2_group.set(gromacs.group_list[gromacs2.group][0])
		v3_field.set(gromacs.force_field_list[gromacs2.field-1][0])
		v4_water.set(gromacs.water_list[gromacs2.water-1][0])
		water_v.set(gromacs.water_list[v4_water.get()-1][1])
		if restraints_var == 0 and config_button_restraints != "":
			checkbox_restraints.set(restraints_var)
			config_button_restraints.configure(state=DISABLED)
		elif restraints_var == 1 and config_button_restraints != "":
			checkbox_restraints.set(restraints_var)
			config_button_restraints.configure(state=ACTIVE)
	else:
		progress = Progress_status(dir_path_project)
		if os.path.isfile(dir_path_dynamics + "em.mdp"):
			shutil.copy(dir_path_dynamics + "em.mdp", dir_path_project + "em.mdp")
			em_file_config = open(dir_path_dynamics + "em.mdp", "r").read()
			em_file = Mdp_config("em.mdp",name, em_file_config, 1)
			print "Found em.mdp file. Using it instead of local configuration."
		else:
			em_file = Mdp_config("em.mdp",name,em_init_config, 0)
		if os.path.isfile(dir_path_dynamics + "pr.mdp"):
			shutil.copy(dir_path_dynamics + "pr.mdp", dir_path_project + "pr.mdp")
			pr_file_config = open(dir_path_dynamics + "pr.mdp", "r").read()
			pr_file = Mdp_config("pr.mdp",name, pr_file_config, 1)
			print "Found pr.mdp file. Using it instead of local configuration."
		else:
			pr_file = Mdp_config("pr.mdp",name,pr_init_config, 0)
		if os.path.isfile(dir_path_dynamics + "md.mdp"):
			shutil.copy(dir_path_dynamics + "md.mdp", dir_path_project + "md.mdp")
			md_file_config = open(dir_path_dynamics + "md.mdp", "r").read()
			md_file = Mdp_config("md.mdp",name, md_file_config, 1)
			print "Found md.mdp file. Using it instead of local configuration."
		else:
			md_file = Mdp_config("md.mdp",name,md_init_config, 0)
	try:
		check_var.set(progress.from_begining)
		if main_bar != 0:
			main_status_bar(check_var.get(), main_bar)
	except:
		pass
	save1()
	if from_pymol == 1:
		print "cmd saved"
		cmd.save(dir_path_project+name+".pdb", name) #PyMOL API

#This function will create the window with configuration files based on MDP class
def mdp_configure(config_name, master):
	
	root2 = Toplevel(master)
	
	if config_name == "em":
		options = em_file.options
	elif config_name == "pr":
		options = pr_file.options
	elif config_name == "md":
		options = md_file.options
	
	if config_name == "em":
		root2.wm_title("Energy Minimalization Options")
	elif config_name == "pr":
		root2.wm_title("Position Restrained MD Options")
	elif config_name == "md":
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

##This function will update MDP class objects alfter closeing "mdp_configure" window
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
##This function will activ or disable restraints button in main window based on check box
def restraints(check, config_button):
	global restraints_var
	if check == 1:
		config_button.configure(state=ACTIVE)
		md_file.options[2][0] = "define"
		md_file.update()
		restraints_var = 1
		gromacs.restraints_index()
	elif check == 0:
		config_button.configure(state=DISABLED)
		md_file.options[2][0] = ";define"
		md_file.update()
		restraints_var = 0

##This function will create restraints window
def restraints_window(master):
	root1 = Toplevel(master)
	root1.wm_title("Restraints Configure")
	
	sw = ScrolledWindow(root1, scrollbar=Y) # just the vertical scrollbar
	sw.pack(fill=BOTH, expand=1)
	
	check_var = IntVar(sw.window)
	check_var.set(gromacs2.restraints_nr)
	
	atom_list = []
	number = 0
	
	for group in gromacs.restraints:
		select = Radiobutton(sw.window, text=group[0], value=number, variable=check_var, command=lambda : modify_index(check_var.get(), atom_list[check_var.get()]))
		select.pack()
		text = Text(sw.window)
		text.insert(END, group[1])
		text.pack()
		atom_list.append(text)
		number = number + 1
	
	if plugin == 1:
		select1 = Radiobutton(sw.window, text="[ PyMol Selected ]", value=number, variable=check_var,command=lambda : modify_index(check_var.get(), atom_list[check_var.get()]))
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
		atom_list.append(text1)
		
	ok_button = Button(sw.window, text="OK", command=lambda : modify_index(check_var.get(), atom_list[check_var.get()], root1))
	ok_button.pack()

##This function will modyfie index_dynamics.ndx file based on user choosed restraints		
def modify_index(index_nr, text, root_to_kill=""):
	gromacs2.restraints_nr = index_nr
	if index_nr < len(gromacs.restraints):
		gromacs.restraints[index_nr][1] = text.get(1.0, END)
	index_file = open("index_dynamics.ndx", "w")
	index_file.write("[ Dynamics Selected ]\n"+text.get(1.0, END))
	index_file.close()
	if root_to_kill != "":
		root_to_kill.destroy()

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
	tkMessageBox.showinfo("Clean", "Temporary files are now removed!")
	clean_option()

##--Comand Line Interface--
def dynamics(name = "h"):
	name = name.lower()
	file_path = dir_path_project + name
	field = str(gromacs2.field) + "\n" + str(gromacs2.water)
	group = str(gromacs2.group)
	global status, stop
	stop = 0
	print "Starting PyMOL plugin 'dynamics' ver."+plugin_ver+" by Tomasz Makarewicz"
	
	##If help is called
	if help_name.count(name) == 1:
		print help_option()
		status = ["fail", "Help printed"]
	##If clean is called
	elif clean_name.count(name) == 1:
		clean_option()
		status = ["fail", "Cleaned"]
	
	##Starting real calculations
	##Preprocess
	if status[0] == "ok":
		preproces(name, file_path)
		save1()
	
	##Checking GROMACS
	if status[0] == "ok":
		status = gromacs.status
	
	##Saving configuration files
	if status[0] == "ok" and stop == 0 and progress.to_do[0] == 1:
		mdp_files(name)
		if status[0] == "ok":
			progress.status[0] = 1
			progress.to_do[0] = 0
			save1()

	##Checking variables - temporary disabled
	#if status[0] == "ok":
	#	status = check_variable(name, field1, field2, group)
	#elif status[0] == "fail":
	#	pass
	
	##Counting topology
	if status[0] == "ok" and stop == 0 and progress.to_do[1] == 1:
		gromacs2.pdb2top(name, file_path, field)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[1] = 1
			progress.to_do[1] = 0
			save1()

	##Adding water box
	if status[0] == "ok" and stop == 0 and progress.to_do[2] == 1:
		gromacs2.waterbox(name, file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[2] = 1
			progress.to_do[2] = 0
			save1()
	
	##EM	
	if status[0] == "ok" and stop == 0 and progress.to_do[3] == 1:
		gromacs2.em(name, file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[3] = 1
			progress.to_do[3] = 0
			save1()
	
	##PR
	if status[0] == "ok" and stop == 0 and progress.to_do[4] == 1:
		gromacs2.pr(name, file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[4] = 1
			progress.to_do[4] = 0
			save1()
	
	##Restraints
	if status[0] == "ok" and stop == 0 and restraints_var == 1 and progress.to_do[5] == 1:
		gromacs2.restraints(name)
		status = gromacs2.status
	
	##MD
	if status[0] == "ok" and stop == 0 and progress.to_do[5] == 1:
		gromacs2.md(name, file_path)
		status = gromacs2.status
		if status[0] == "ok":
			progress.status[5] = 1
			progress.to_do[5] = 0
			save1()
	
	##Trjconv
	if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
		gromacs2.trjconv(name, file_path, group)
		status = gromacs2.status
	
	##Showing multimodel
	if status[0] == "ok" and stop == 0 and progress.to_do[6] == 1:
		show_multipdb(name)
		progress.status[6] = 1
		progress.to_do[6] = 0
		save1()
	elif status[0] == "fail":
		print status[1]
		if help_name.count(name) != 1 and clean_name.count(name) != 1:
			error_message()
		status = ["ok", ""]
	#This should be removed
	if plugin == 0:
		return dir_path_project

##Preprocesing
def preproces(name, file_path):
	status = ["ok", ""]
	##Checking directories
	if os.path.isdir(dir_path_project) == False:
		os.makedirs(dir_path_project)
	os.chdir(dir_path_project)

##Checking if given varaibles are correct - depreciated
def check_variable(name, field1, field2, group):
	#jeszcze nie dziala (name)
	status = ["ok", ""]
	try:
		int(field1)
	except:
		status = ["fail", "Wrong Force Field"]
	try:
		int(field2)
	except:
		status = ["fail", "Wrong Water Model"]
	try:
		int(group)
	except:
		status = ["fail", "Wrong Group"]
	return status

##Saving configuration files
def mdp_files(name):
	em_file.save_file(dir_path_project)
	pr_file.save_file(dir_path_project)
	md_file.save_file(dir_path_project)

##Show multimodel PDB file in PyMOOL
def show_multipdb(name):
	if plugin == 1:
		try:
			cmd.hide("everything", name) #PyMOL API
		except:
			pass
		cmd.load(name+"_multimodel.pdb") #PyMOL API

##Saving tar.bz file
def save(destination_path):
	print "Saving"
	import tarfile
	save1()
	name = dir_path_project.split("/")
	preproces(name[-2], dir_path_project+name[-2])
	tar = tarfile.open(destination_path+".tar.bz2", "w:bz2")
	tar.add(dir_path_project, recursive=True, arcname=name[-2])
	tar.close()
	os.remove(destination_path)

##Load tar.bz file
def load(file_path):
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
	global dir_path_project
	dir_path_project = dir_path_dynamics + names[0] + "/"
	load1()
	return names[0]

##Save all settings to options.pickle file
def save1():
	if os.path.isdir(dir_path_project) == False:
		os.makedirs(dir_path_project)
	destination_option = file(dir_path_project + "options.pickle", "w")
	pickle.dump([plugin_ver, gromacs.version, gromacs2, em_file, pr_file, md_file, progress, restraints_var], destination_option)
	del destination_option
	#print "options saved"

##Load all settings from options.pickle file	
def load1():
	global gromacs, gromacs2, em_file, pr_file, md_file, progress, restraints_var
	
	pickle_file = file(dir_path_project +"options.pickle")
	options = pickle.load(pickle_file)
	
	print "Project was created for Dynamics PyMOL Plugin"+options[0]+" and "+options[1]
	if gromacs.version != options[1]:
		print "GROMACS versions not identical."
	if options[0][1:4] == "1.0":
		print "1.0 compatibility layer"
		#gromacs.version = options[1]
		gromacs2 = options[2]
		em_file = options[3]
		pr_file = options[4]
		md_file = options[5]
		progress = options[6]
		restraints_var = options[7]
		if restraints_var == 1:
			gromacs.restraints_index()

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
Choose force field and water model options in the middle column.
Choose any additional options in the right column.
Press OK button.
Click Start button and wait till calculations are finished.
You will see animated model in PyMOL viewer."""
	return help_message

##Clean function
def clean_option():
	dir_path = os.getenv("HOME")+"/.dynamics/"
	shutil.rmtree(dir_path)
	print "Temporary files are now removed!"

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
	
##--Comand Line Interface-- - depreciated
if plugin == 1:
	cmd.extend("dynamics", dynamics)
