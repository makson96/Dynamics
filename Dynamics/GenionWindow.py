import pymol_plugin_dynamics
from tkinter import *
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
        ok_button = Button(root, text="OK", command=lambda: [gromacs2.update(
            {"salt_conc": salt.get(), "positive_ion": posit.get(), "negative_ion": negat.get()}), root.destroy()])
        ok_button.pack(side=TOP)
