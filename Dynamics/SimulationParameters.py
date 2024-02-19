import GromacsOutput
import GromacsInput

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
            self.vectors_prody = Vectors()
            print("ProDy correctly imported")
        self.progress = ProgressStatus()

    def create_cfg_files(self):
        self.em_file, self.pr_file, self.md_file = create_config_files(self.project_name)

    def change_stop_value(self, value):
        if value:
            self.stop = True
        else:
            self.stop = False

    def change_project_name(self, name):
        self.project_name = name
        project_dir = get_project_dirs(self.project_name)
        if not os.path.isdir(project_dir):
            os.makedirs(project_dir)

