import pymol_plugin_dynamics
import queue as Queue
from tkinter import *
import _thread as thread
from pymol import cmd, cgo, parsing, plugins, CmdException
try:
    import prody
except ModuleNotFoundError:
    prody = False
import time
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
        project_dir = pymol_plugin_dynamics.get_project_dirs(project_name)
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
        vectors_prody = s_params.vectors_prody
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
        save_button = Button(frame1_8, text="Save", command=lambda: pymol_plugin_dynamics.select_file_save(s_params))
        save_button.pack(side=LEFT)
        log_button = Button(frame1_8, text="Log", command=pymol_plugin_dynamics.log_window(s_params))
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
