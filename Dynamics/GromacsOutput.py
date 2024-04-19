# This class is responsible for interface to GROMACS. It will read all important data from GROMACS tools.
import pymol_plugin_dynamics
import os
import re


class GromacsOutput:

    def __init__(self):
        # Remove garbage
        self.dynamics_dir = pymol_plugin_dynamics.get_dynamics_dir()
        self.gmx_exe = ""
        self.version = ""
        self.build_arch = ""
        self.build_on_cygwin = 0
        self.get_gromacs_exe_info()

        self.force_list = []
        self.water_list = []
        self.group_list = []
        self.restraints = []
        if self.gmx_exe:
            # Track current directory and switch to dynamics_dir before invoking gmx...
            current_dir = os.getcwd()
            os.chdir(self.dynamics_dir)
            self.gmx_stdout_file_path = "test_gromacs.txt"
            self.create_gmx_stdout_file()
            self.init2()
            # Switch back to current directory...
            os.chdir(current_dir)

    def create_gmx_stdout_file(self):
        garbage_files = next(os.walk(self.dynamics_dir))[2]
        for garbage in garbage_files:
            if garbage[0] == "#":
                os.remove(self.dynamics_dir + garbage)

        print("Reading available force fields and water models")
        tg_pdb_file_path = "test_gromacs.pdb"
        tg_pdb_file = open(tg_pdb_file_path, "wb")
        tg_pdb_file.write(b"ATOM      1  N   LYS     1      24.966  -0.646  22.314  1.00 32.74      1SRN  99\n")
        tg_pdb_file.close()

        gmx_stdin_file_path = "gromacs_stdin.txt"
        gmx_stdin_file = open(gmx_stdin_file_path, "w")
        gmx_stdin_file.write("1\n1")
        gmx_stdin_file.close()

        cmd = f"{self.gmx_exe} pdb2gmx -f {tg_pdb_file_path} -o test_gromacs.gro -p test_gromacs.top"
        pymol_plugin_dynamics.execute_subprocess(cmd, gmx_stdin_file_path, self.gmx_stdout_file_path)

    def init2(self):
        lista_gromacs = pymol_plugin_dynamics.read_lines_standardize(self.gmx_stdout_file_path)

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
        self.water_list = pymol_plugin_dynamics.get_water_models_info(lista_gromacs)

        print("Reading available groups")

        gmx_stdin_file_path = "gromacs_stdin.txt"
        fo = open(gmx_stdin_file_path, "w")
        fo.write("1")
        fo.close()

        gmx_stdout_file_path = "test_gromacs.txt"

        cmd = "{} trjconv -f  test_gromacs.pdb -s test_gromacs.pdb -o test_gromacs2.pdb".format(self.gmx_exe)
        pymol_plugin_dynamics.execute_subprocess(cmd, gmx_stdin_file_path, gmx_stdout_file_path)

        group_test_list = pymol_plugin_dynamics.read_lines_standardize(gmx_stdout_file_path)

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
        os.chdir(pymol_plugin_dynamics.get_dynamics_dir())

        print("Updating available water models")
        gmx_stdin_file_path = "gromacs_stdin.txt"
        fo = open(gmx_stdin_file_path, "w")
        fo.write("%d\n" % force_number)
        fo.write("1")
        fo.close()

        gmx_stdout_file_path = "test_gromacs.txt"

        cmd = "{} pdb2gmx -f  test_gromacs.pdb -o test_gromacs.gro -p test_gromacs.top".format(self.gmx_exe)
        pymol_plugin_dynamics.execute_subprocess(cmd, gmx_stdin_file_path, gmx_stdout_file_path)

        lista_gromacs = pymol_plugin_dynamics.read_lines_standardize(gmx_stdout_file_path)
        self.water_list = pymol_plugin_dynamics.get_water_models_info(lista_gromacs)

        # Switch back to current directory...
        os.chdir(current_dir)

        # save_options()
        return self.water_list

    # This function will read atoms group for restraints for current molecule.
    def restraints_index(self, project_name):
        self.restraints = []
        current_dir = os.getcwd()
        os.chdir(pymol_plugin_dynamics.get_project_dirs(project_name))

        fo = open("gromacs_stdin.txt", "w")
        fo.write("q")
        fo.close()

        cmd = "{} make_ndx -f {}.pdb -o index.ndx".format(self.gmx_exe, project_name)
        pymol_plugin_dynamics.execute_subprocess(cmd, "gromacs_stdin.txt", "restraints.log")

        index_list = pymol_plugin_dynamics.read_lines_standardize("restraints.log")

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

    # Detect gmx executable along with other associated information...
    def get_gromacs_exe_info(self):
        gmx_exes = ['gmx_mpi_d', 'gmx_mpi', 'gmx']

        for gmx in gmx_exes:
            cmd = f"{gmx} -version"
            output = pymol_plugin_dynamics.check_output_subprocess(cmd)
            if not re.search("GROMACS version:", output, re.I):
                continue

            self.gmx_exe = gmx
            for line in output.split("\n"):
                if re.search("^[ ]*GROMACS version:", line, re.I):
                    self.gmx_exe = gmx
                    self.version = re.sub("^[ ]*GROMACS version:[ ]*", "", line, flags=re.I)
                    if "VERSION " in self.version:
                        self.version = self.version.split("VERSION ")[1].rstrip()
                elif re.search(r"^[ ]*Build OS/arch:", line, re.I):
                    self.build_arch = re.sub("^[ ]*Build OS/arch:[ ]*", "", line, flags=re.I)

                    if re.search(r"CYGWIN", self.build_arch, re.I):
                        self.build_on_cygwin = 1
            break
