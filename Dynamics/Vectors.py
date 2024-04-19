import pymol_plugin_dynamics
from tkinter import *
try:
    import prody
except ModuleNotFoundError:
    prody = False
from pymol import cmd, cgo, parsing, plugins, CmdException
# This class will handle PCA by ProDy python library and show vectors from NMD file.
class Vectors:
    nmd_name = []
    nmd_atomnames = []
    nmd_resnames = []
    nmd_resids = []
    nmd_bfactors = []
    nmd_coordinates = []
    nmd_mode = []
    nmd_scale_mode = []
    color = "grey"
    scale = 1.0
    mode_nr = 0

    calculation_type = 0
    contact_map = 0
    block_contact_map = 0

    enm = 0

    # Change Multimodel PDB file into NMD vector file
    def prody(self, project_name):
        # Silence ProDy and create logs
        prody.confProDy(verbosity='none')
        prody.startLogfile("log_prody.log")
        # Prepare ensemble
        model = prody.parsePDB(project_name + "_multimodel.pdb", subset='calpha')
        ensemble = prody.Ensemble(project_name + ' ensemble')
        ensemble.setCoords(model.getCoords())
        ensemble.addCoordset(model.getCoordsets())
        ensemble.iterpose()
        # ANM calculations
        if self.calculation_type == 0:
            anm = prody.ANM(project_name)
            anm.buildHessian(ensemble)
            anm.calcModes()
            write_nmd = anm
            self.enm = anm
        # PCA calculations
        elif self.calculation_type == 1:
            pca = prody.PCA(project_name)
            pca.buildCovariance(ensemble)
            pca.calcModes()
            write_nmd = pca
        # GNM calculations
        elif self.calculation_type == 2:
            gnm = prody.GNM(project_name)
            gnm.buildKirchhoff(ensemble)
            gnm.calcModes()
            write_nmd = gnm
            self.enm = gnm
        # Write NMD file
        prody.writeNMD(project_name + '.nmd', write_nmd[:3], model)
        prody.closeLogfile("log_prody.log")

    # Read NMD file
    def nmd_format(self, project_name):
        file_nmd = open('{}.nmd'.format(project_name), "r")
        list_nmd = file_nmd.readlines()

        self.nmd_mode = []
        self.nmd_scale_mode = []
        for line in list_nmd:
            split_line = line.split()
            if split_line[0] == "name":
                self.nmd_name = split_line
                self.nmd_name.pop(0)
            elif split_line[0] == "atomnames":
                self.nmd_atomnames = split_line
                self.nmd_atomnames.pop(0)
            elif split_line[0] == "resnames":
                self.nmd_resnames = split_line
                self.nmd_resnames.pop(0)
            elif split_line[0] == "resids":
                self.nmd_resids = split_line
                self.nmd_resids.pop(0)
            elif split_line[0] == "bfactors":
                self.nmd_bfactors = split_line
                self.nmd_bfactors.pop(0)
            elif split_line[0] == "coordinates":
                self.nmd_coordinates = split_line
                self.nmd_coordinates.pop(0)
            elif split_line[0] == "mode":
                pre_mode = split_line
                self.nmd_mode.append(pre_mode[3:])
                self.nmd_scale_mode.append(pre_mode[2])

    # Show contact map on PyMOL screen
    def show_contact_map(self, sensitivity, project_name):
        contact_matrix = self.enm.getKirchhoff()
        c_alpha_nr = 0
        for c_alpha_list in contact_matrix:
            c_alpha_nr = c_alpha_nr + 1
            c_alpha_target_nr = 0
            for c_alpha_1 in c_alpha_list:
                c_alpha_target_nr = c_alpha_target_nr + 1
                if c_alpha_nr != c_alpha_target_nr and float(c_alpha_1) < float(sensitivity):
                    cmd.select("sele1",
                               "n. ca and {}_multimodel and i. {}".format(project_name, str(c_alpha_nr)))  # PyMOL API
                    cmd.select("sele2", "n. ca and {}_multimodel and i. {}".format(project_name,
                                                                                   str(c_alpha_target_nr)))  # PyMOL API
                    cmd.distance("contact_map", "sele1", "sele2")  # PyMOL API
        try:
            cmd.hide("labels", "contact_map")  # PyMOL API
            cmd.delete("sele1")  # PyMOL API
            cmd.delete("sele2")  # PyMOL API
        except:
            pass

    # Show contact map/cross corelation as a graph
    def graph_contact_map(self, plot_type):
        if plot_type == "contact":
            # matplotlib
            prody.showContactMap(self.enm)
        elif plot_type == "cross":
            # matplotlib
            prody.showCrossCorr(self.enm)

    # Show vectors from NMD file
    def show_vectors(self):
        color1 = cmd.get_color_tuple(self.color)  # PyMOL API
        color2 = cmd.get_color_tuple(self.color)  # PyMOL API
        if color1:
            color1 = list(color1)
        # Fallback to grey in case of unrecognized color
        else:
            color1 = [0.5, 0.5, 0.5]
        if color2:
            color2 = list(color2)
        # Fallback to grey in case of unrecognized color
        else:
            color2 = [0.5, 0.5, 0.5]
        arrow_head_radius = 0.15

        x1 = []
        y1 = []
        z1 = []

        coor = "x"
        for coordinate in self.nmd_coordinates:
            if coor == "x":
                x1.append(float(coordinate))
                coor = "y"
            elif coor == "y":
                y1.append(float(coordinate))
                coor = "z"
            elif coor == "z":
                z1.append(float(coordinate))
                coor = "x"

        x2 = []
        y2 = []
        z2 = []

        # This factor is provided to make vector length more like in NMWiz.
        # More investigation is needed to get exact formula.
        approximation_factor = 16.6

        coor = "x"
        coor_nr = 0
        round_nr = 0
        for mode in self.nmd_mode[self.mode_nr]:
            if coor == "x":
                x2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + x1[
                        coor_nr])
                coor = "y"
            elif coor == "y":
                y2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + y1[
                        coor_nr])
                coor = "z"
            elif coor == "z":
                z2.append(
                    float(mode) * float(self.nmd_scale_mode[self.mode_nr]) * approximation_factor * self.scale + z1[
                        coor_nr])
                coor = "x"
            round_nr = round_nr + 1
            if round_nr == 3:
                round_nr = 0
                coor_nr = coor_nr + 1

        coor_nr = 0
        for position in x1:
            try:
                cmd.delete("Mode_Vector_" + str(coor_nr))
            except:
                pass

            cone = [cgo.CONE, x1[coor_nr], y1[coor_nr], z1[coor_nr], x2[coor_nr], y2[coor_nr], z2[coor_nr],
                    arrow_head_radius, 0.0] + color1 + color2 + [1.0, 0.0]
            cmd.load_cgo(cone, "Mode_Vector_" + str(coor_nr))  # PyMOL API
            coor_nr = coor_nr + 1
        # Another workaround for PyMOL 1.8 with TravisCI
        try:
            cam_possition = cmd.get_view(quiet=1)  # PyMOL API
            cmd.set_view(cam_possition)  # PyMOL API
        except TypeError:
            pass

    def change_vectors_color(self, color):
        self.color = color
        self.show_vectors()

    def change_vectors_scale(self, scale):
        scale = float(scale)
        self.scale = scale
        self.show_vectors()

    def change_vectors_mode_nr(self, mode_nr):
        self.mode_nr = mode_nr
        self.show_vectors()

    def options_change(self, v1, v2, root):
        self.calculation_type = v1.get()
        self.contact_map = v2.get()
        # save_options()
        root.destroy()

    def block_contact(self, block, contact_map_b, contact_map_v):
        # ToDO: Replace below Tk code with something agnostic
        self.block_contact_map = block
        if block == 0:
            contact_map_b.configure(state=ACTIVE)
        elif block == 1:
            contact_map_b.configure(state=DISABLED)
            contact_map_v.set(0)

