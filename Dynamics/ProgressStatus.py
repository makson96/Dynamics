import pymol_plugin_dynamics
# Status and to_do maintaining class
class ProgressStatus:
    # 0:Save configuration files; 1:Generate topology file from pdb; 2:Adding Water Box;
    # 3: Adding ions and neutralization 4:Energy Minimization; 5:Position Restrained MD;
    # 6:Restraints; 7:Molecular Dynamics Simulation; 8:Generate multimodel PDB; 9:Calculate vectors using ProDy
    status = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    to_do = [1, 1, 1, 1, 1, 1, 0, 1, 1, 1]

    resume = 0
    x2top = 0
    steps = 8

    def to_do_update(self, position, value):
        if isinstance(position, int):
            self.to_do[position] = value
            self.to_do = self.to_do

    #            save_options()

    def x2top_update(self, value):
        if isinstance(value, int):
            self.x2top = value

    def to_do_status(self):
        to_do = []
        for work in self.status:
            if work == 0:
                to_do.append(1)
            elif work == 1:
                to_do.append(0)
        self.to_do = to_do
