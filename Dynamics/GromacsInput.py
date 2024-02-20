import pymol_plugin_dynamics
import os
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
        pymol_plugin_dynamics.status_update(status)
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
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        if os.path.isfile("{}.gro".format(project_name)):
            status = ["ok", ""]
        else:
            status = ["fail", "Warning. Trying to ignore unnecessary hydrogen atoms."]

            command = "{0} pdb2gmx -ignh -f {1}.pdb -o {1}.gro -p {1}.top {2}".format(gmx_cmd, project_name, hh)
            pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')
        pymol_plugin_dynamics.status_update(status)

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
        pymol_plugin_dynamics.status_update(status)
        gmx_cmd = s_params.gmx_output.command
        project_name = s_params.project_name
        try:
            os.remove("{}.gro".format(project_name))
            os.remove("{}.top".format(project_name))
        except FileNotFoundError:
            pass

        command = "{0} x2top -f {1}.pdb -o {1}.top".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}.top".format(project_name)) and not stop:
            status = ["ok", "Calculating structure using trjconv."]
        else:
            status = ["fail", "Unable to create topology file."]
        pymol_plugin_dynamics.status_update(status)

        if status[0] == "ok":
            fo = open("gromacs_stdin.txt", "w")
            fo.write("0")
            fo.close()

            command = "{0} trjconv -f {1}.pdb -s {1}.pdb -o {1}.gro".format(gmx_cmd, project_name)
            pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} editconf -f {1}.gro -o {1}1.gro -c {2} {3} {4}".format(gmx_cmd, project_name, box_type, distance,
                                                                              density)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

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
        pymol_plugin_dynamics.status_update(status)

        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

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
        pymol_plugin_dynamics.status_update(status)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        fo = open("gromacs_stdin.txt", "w")
        fo.write("13")
        fo.close()

        status = ["ok", "Adding salts and ions"]
        pymol_plugin_dynamics.status_update(status)

        command = "{0} genion -s {1}_ions.tpr -o {1}_b4em.gro {2} {3} {4} {5} -p {1}.top".format(gmx_cmd,
                                                                                                 project_name, positive,
                                                                                                 negative, salt, neu)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} grompp -f em -c {1}_b4em -p {1} -o {1}_em".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_em -o {1}_em -c {1}_b4pr -v".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} grompp -f pr -c {1}_b4pr -r {1}_b4pr -p {1} -o {1}_pr".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_pr -o {1}_pr -c {1}_b4md -v".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} genrestr -f {1}.pdb -o posre_2.itp -n index_dynamics.ndx".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} grompp -f md -c {1}_b4md  -p {1} -o {1}_md".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

        command = "{0} mdrun -nice 4 -s {1}_md -o {1}_md -c {1}_after_md -v".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, None, 'log1.txt', 'log.txt')

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

        pymol_plugin_dynamics.status_update(status)
        command = "{0} trjconv -f {1}_md.trr -s {1}_md.tpr -o {1}_multimodel.pdb".format(gmx_cmd, project_name)
        pymol_plugin_dynamics.execute_and_monitor_subprocess(command, 'gromacs_stdin.txt', 'log1.txt', 'log.txt')

        stop = s_params.stop
        if os.path.isfile("{}_multimodel.pdb".format(project_name)) and not stop:
            status = ["ok", "Finished!"]
        else:
            status = ["fail", "Unable to generate multimodel PDB file"]
        return status
