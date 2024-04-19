# This class create and maintain abstraction mdp file representatives. em.mdp, pr.mdp, md.mdp

import pymol_plugin_dynamics
class MdpConfig:
    external_file = 0
    options = [[]]
    file_name = ""

    def __init__(self, file_name, init_config, external_file=0):
        self.file_name = file_name
        self.external_file = external_file
        list1 = init_config.split("\n")
        list2 = []
        for line in list1:
            list2.append(line.split(" = "))
        self.options = list2

    def update(self, option_nr, value, check=1):
        self.options[option_nr][1] = value
        if check == 0 and self.options[option_nr][0][0] != ";":
            self.options[option_nr][0] = ";" + self.options[option_nr][0]
        elif check == 1 and self.options[option_nr][0][0] == ";":
            self.options[option_nr][0] = self.options[option_nr][0][1:]
        self.clean_artefacts()

    def save_file(self, s_params):
        project_name = s_params.project_name
        project_dir = pymol_plugin_dynamics.get_project_dirs(project_name)
        config = ""
        for option in self.options:
            # pass empty option
            if option == ['']:
                pass
            else:
                config = "{}{} = {}\n".format(config, str(option[0]), str(option[1]))
        mdp = open(project_dir + self.file_name, "w")
        mdp.write(config)
        mdp.close()

    # Clean options from artefacts
    def clean_artefacts(self):
        try:
            self.options.remove([''])
        except:
            pass
