import GromacsOutput
import GromacsInput
import pymol_plugin_dynamics
import ProgressStatus
import os
import Vectors
try:
    import prody
except ModuleNotFoundError:
    prody = False
class SimulationParameters:
    gmx_output = ""
    gmx_input = ""
    vectors_prody = False
    stop = False
    project_name = "nothing"
    progress = ""
    em_file = ""
    pr_file = ""
    md_file = ""

    def __init__(self):
        self.gmx_output = GromacsOutput.GromacsOutput()
        self.gmx_input = GromacsInput.GromacsInput()
        print("Found GROMACS VERSION {}".format(self.gmx_output.version))
        if prody:
            self.vectors_prody = Vectors.Vectors()
            print("ProDy correctly imported")
        self.progress = ProgressStatus.ProgressStatus()

    def create_cfg_files(self):
        self.em_file, self.pr_file, self.md_file = pymol_plugin_dynamics.create_config_files(self.project_name)

    def change_stop_value(self, value):
        if value:
            self.stop = True
        else:
            self.stop = False

    def change_project_name(self, name):
        self.project_name = name
        project_dir = pymol_plugin_dynamics.get_project_dirs(self.project_name)
        if not os.path.isdir(project_dir):
            os.makedirs(project_dir)

