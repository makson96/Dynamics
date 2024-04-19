import pymol_plugin_dynamics
import os
from tkinter import *
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
        self.change_e(gromacs2.explicit, v4_water, water_v, v2, s_params)

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
        dynamics_dir = pymol_plugin_dynamics.get_dynamics_dir()
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
        dynamics_dir = pymol_plugin_dynamics.get_dynamics_dir()
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
            {"box_distance": distance.get(), "box_density": density.get()}))
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
