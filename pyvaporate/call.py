import subprocess

import os


LAMMPS_CMD = os.environ["LAMMPS_CMD"]
TAPSIM_CMD = os.environ["TAPSIM_BIN"] + "/tapsim"
MESHGEN_CMD = os.environ["TAPSIM_BIN"] + "/meshgen"
MGN_INI_LINES = open("{}/meshgen.ini".format(
    "/".join(__file__.split("/")[:-1]))
).readlines()
ELTS = {"10": "W"}
E_FIELDS = {"10": "57.1e-9"}
MASSES = {"W": "183.85"}
CHARGE_STATES = {"W": "3"}

def call_tapsim(node_file, n_events):
    """
    Run the TAPSim program, starting with building a voronoi mesh for
    an emitter (node) file and then running a set number of
    evaporation steps before stopping.
    """

    print("Writing meshgen.ini")
    write_meshgen_ini()

    print("Creating sampleMesh.bin and sampleMesh.cfg")
    _ = subprocess.check_output(
        [MESHGEN_CMD, node_file, "sampleMesh.bin",
         "--create-config-template=sampleMesh.cfg"]
    )

    sm_lines = open("sampleMesh.cfg").readlines()
    with open("sampleMesh.cfg", "w") as sm:
        edit = False
        for line in sm_lines:
            if "ID" in line and line.split()[-1] not in ["0", "1", "2", "3"]:
                ID = line.split()[-1]
                elt = ELTS[ID]
                edit = True
            if edit:
                if "***_SET_NAME_HERE_***" in line:
                    line = line.replace("***_SET_NAME_HERE_***", elt)
                elif "***_SET_MASS_HERE_***" in line:
                    line = line.replace("***_SET_MASS_HERE_***", MASSES[elt])
                elif "EVAPORATION_CHARGE_STATE" in line:
                    line = line.replace("1", CHARGE_STATES[elt])
                elif "***_SET_EVAPORATION_FIELD_STRENGTH_HERE_***" in line:
                    line = line.replace(
                        "***_SET_EVAPORATION_FIELD_STRENGTH_HERE_***",
                        E_FIELDS[ID]
                    )
                sm.write(line)
            else:
                sm.write(line)

    print("Running TAPSim")
    _ = subprocess.check_output(
        [TAPSIM_CMD, "evaporation", "sampleMesh.cfg", "sampleMesh.bin",
         "--event-limit={}".format(n_events), "--write-ascii"]
    )


def update_emitter(node_file, results_file):
    """
    """
    emitter_lines = open(node_file).readlines()
    results_lines = open(results_file).readlines()
    results_lines = results_lines[results_lines.index("ASCII\n")+1:]
    remove_numbers = [line.split()[2] for line in results_lines]
    emitter_lines = [
        l for l in emitter_lines if l.split()[-1] not in remove_numbers
    ]
    with open(node_file, "w") as nf:
        for line in emitter_lines:
            nf.write(line)

def write_meshgen_ini():
    with open("meshgen.ini", "w") as mgn:
        for line in MGN_INI_LINES:
            mgn.write(line)


def relax_emitter(emitter_file):

    #TODO: Build emitter_relax.in
    convert_emitter_to_lammps(emitter_file)
    write_lammps_input_file("data.emitter")

    _ = subprocess.check_output([LAMMPS_CMD, "in.emitter_relax"])


def find_surface_atoms():
    surface_file = [f for f in os.listdir(os.getcwd()) if "surface_data" in f][-1]
    surface_lines = open(surface_file).readlines()
    surface_numbers = []
    for line in surface_lines[5:]:
        split_line = line.split()
        if split_line[1] == "10":
            surface_numbers.append(split_line[0])
    return surface_numbers


def convert_emitter_to_lammps(emitter_file, surface_numbers):

    emitter_lines = open(emitter_file).readlines()
    atom_lines = [l for l in emitter_lines[1:-1] if l.split()[-2] not in
                  ["0", "1", "2", "3"]]
    atom_types = [t.split("=")[-1] for t in emitter_lines[-1].split()[1:]]
    print(atom_types)
    atom_coords = []
    for line in atom_lines:
        split_line = line.split()
        atom_coords.append(
            [split_line[-1], str(atom_types.index(ELTS[split_line[-2]])+1)]+
            [str(float(x)*1e10) for x in split_line[:3]]
        )

    x_coords = [float(a[1]) for a in atom_coords]
    y_coords = [float(a[2]) for a in atom_coords]
    z_coords = [float(a[3]) for a in atom_coords]
    xlim = (min(x_coords), max(x_coords))
    ylim = (min(y_coords), max(y_coords))
    zlim = (min(z_coords), max(z_coords))

    with open("data.emitter", "w") as dat:
        dat.write("LAMMPS Emitter\n\n")
        dat.write("{} atoms\n\n".format(len(atom_lines)))
        dat.write("{} atom types\n\n".format(len(atom_types)))
        dat.write("{} {} xlo xhi\n".format(xlim[0], xlim[1]))
        dat.write("{} {} ylo yhi\n".format(ylim[0], ylim[1]))
        dat.write("{} {} zlo zhi\n\n".format(zlim[0], zlim[1]))
        dat.write("Masses\n\n")
        for elt in atom_types:
            dat.write("{} {}\n".format(atom_types.index(elt)+1, MASSES[elt]))
        dat.write("\nAtoms\n\n")
        i = 1
        surface_indices = []
        for atom in atom_coords:
            if int(atom[0]) in surface_numbers:
                surface_indices.append(i)
            dat.write("{}\n".format(" ".join([str(i)]+atom[1:])))
            i += 1
    with open("surface_indices.txt", "w") as si:
        for index in surface_indices:
            si.write("{}\n".format(index))

def write_lammps_input_file(structure_file):
    with open("in.emitter_relax", "w") as er:
        er.write("# Emitter Relaxation\n\n")
        er.write("units real\n atom_style atomic\n\nread_data {}\n\n".format(structure_file))
        er.write("boundary s s s\n\n")
        er.write("pair_style meam\n")
        er.write("pair_coeff * * /u/mashton/software/lammps/library.meam W NULL W\n\n")
        er.write("neighbor 1.0 bin\n")
        er.write("neigh_modify delay 5 every 1\n\n")
        er.write("fix 1 all nve\n")
        er.write("timestep 0.005\n")
        er.write("run 100")
