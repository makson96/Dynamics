import pymol_plugin_dynamics
import queue as Queue
from tkinter import *
from tkinter.ttk import Progressbar, Scrollbar
import _thread as thread
import time
from tkinter import messagebox as tkMessageBox
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
            pymol_plugin_dynamics.no_molecule_warning()

    # This function will create main Calculation Window
    def window(self, root, s_params, status, parent):
        root.wm_title("Calculation Window")
        frame1 = Frame(root)
        frame1.pack(side=TOP)
        frame2 = Frame(root)
        frame2.pack(side=TOP)
        s_params.change_stop_value(True);
        self.bar_var = StringVar(root)
        self.bar_var.set("Ready to start")

        w5 = Label(frame1, textvariable=self.bar_var)
        w5.pack(side=TOP)
        self.bar_widget = Progressbar(frame1)
        self.bar_widget.pack(side=TOP)

        exit_button = Button(frame2, text="EXIT", command=root.destroy)
        exit_button.pack(side=LEFT)

        save_button = Button(frame2, text="SAVE", command=lambda: pymol_plugin_dynamics.select_file_save(1))
        save_button.pack(side=LEFT)

        stop_button = Button(frame2, text="STOP", command=lambda: self.start_counting(0, s_params))
        stop_button.pack(side=LEFT)
        stop = s_params.stop
        if stop:
            stop_button.configure(state=DISABLED)
        self.stop_button = stop_button

        start_button = Button(frame2, text="START", command=lambda: self.start_counting(1, s_params))
        start_button.pack(side=LEFT)
        if stop == 0:
            start_button.configure(state=DISABLED)
        self.start_button = start_button

        log_button = Button(frame2, text="LOG", command=pymol_plugin_dynamics.log_window)
        log_button.pack(side=LEFT)
        log_button.configure(state=DISABLED)
        self.log_button = log_button

        # Updateing status bar
        tasks_nr = 0.0
        for task in s_params.progress.to_do:
            tasks_nr = tasks_nr + task
        self.tasks_to_do = tasks_nr       
        thread.start_new_thread(self.bar_update, (s_params, status))
        self.bar_display(root, parent, s_params)

    # This function will update status bar during molecular dynamics simulation (beware this is separate thread)
    def bar_update(self, s_params, status):
        percent = 0.0
        while s_params.stop:
            time.sleep(0.5)
        while percent != 100:  # and error == ""
            time.sleep(0.5)
            percent = pymol_plugin_dynamics.steps_status_bar("only_bar", s_params)
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
            self.start_counting(0, s_params)
            self.start_button.configure(state=DISABLED)
            tkMessageBox.showerror("GROMACS Error Message", "Error")  # error)
        if status == "Finished!":
            root.destroy()
            # Show interpretation window after successful completion of the calculations...
            pymol_plugin_dynamics.show_interpretation_window(parent, s_params)
        else:
            root.after(100, self.bar_display, root, parent, s_params)

    # This function will change global value if stop is clicked during simulation
    def start_counting(self, value, s_params):
        if value == 1:
            stop = 0
            thread.start_new_thread(pymol_plugin_dynamics.dynamics, (s_params,))
            self.stop_button.configure(state=ACTIVE)
            self.start_button.configure(state=DISABLED)
            self.log_button.configure(state=DISABLED)
        elif value == 0:
            stop = 1
            self.stop_button.configure(state=DISABLED)
            self.start_button.configure(state=ACTIVE)
            self.log_button.configure(state=ACTIVE)
