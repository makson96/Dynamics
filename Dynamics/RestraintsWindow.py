import pymol_plugin_dynamics
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
