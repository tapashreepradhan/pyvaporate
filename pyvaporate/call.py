import subprocess

import os

from pyvaporate.evaluate import assign_ids_by_cn

from monty.serialization import loadfn


CONFIGURATION = loadfn("/u/mashton/pv.yaml")
LAMMPS_CMD = CONFIGURATION["LAMMPS_CMD"]
TAPSIM_CMD = CONFIGURATION["TAPSIM_CMD"]
MESHGEN_CMD = CONFIGURATION["MESHGEN_CMD"]
MGN_INI_LINES = open("{}/meshgen.ini".format(
    "/".join(__file__.split("/")[:-1]))
).readlines()
ELTS = {"10": "W"}
E_FIELDS = {"10": "57.1e-9", "11": "27.1e-9", "12": "37.1e-9", "13": "47.1e-9", "14": "57.1e-9", "15": "67.1e-9", "16": "77.1e-9", "17": "87.1e-9", "18": "97.1e-9"}
MASSES = {"1": "183.85"}
CHARGE_STATES = {"W": "3"}

def call_tapsim(node_file, n_events, id_dict):
    """
    Run the TAPSim program, starting with building a voronoi mesh for
    an emitter (node) file and then running a set number of
    evaporation steps before stopping.
    """

    write_meshgen_ini()

    print("Creating mesh.txt and mesh.cfg")
    _ = subprocess.check_output(
        [MESHGEN_CMD, node_file, "mesh.txt",
         "--create-config-template=mesh.cfg", "--write-ascii"]
    )

    assign_labels_and_unique_ids(id_dict)

    sm_lines = open("mesh.cfg").readlines()
    with open("mesh.cfg", "w") as sm:
        edit = False
        for line in sm_lines:
            if "ID" in line and line.split()[-1] not in ["0", "1", "2", "3"]:
                ID = line.split()[-1]
                BASE_ID = str(int(ID) - int(ID[-1]))
                elt = ELTS[BASE_ID]
                name = elt+"_"+ID[-1]
                edit = True
            if edit:
                if "***_SET_NAME_HERE_***" in line:
                    line = line.replace("***_SET_NAME_HERE_***", name)
                elif "***_SET_MASS_HERE_***" in line:
                    line = line.replace("***_SET_MASS_HERE_***", MASSES[BASE_ID])
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
        [TAPSIM_CMD, "evaporation", "mesh.cfg", "mesh.txt",
         "--event-limit={}".format(n_events), "--write-ascii"]
    )


def assign_labels_and_unique_ids(id_dict):
    m_lines = open("mesh.txt").readlines()
    i = 1
    with open("mesh.txt", "w") as m:
        m.write(m_lines[0])
        for line in m_lines[1:]:
            sl = line.split()
#            sl[4] = str(i)
            m.write("	".join(sl))
            m.write("\n")
            i += 1
        comment = "#"
        for ID in id_dict:
            comment += " {}={}".format(ID, id_dict[ID])
        m.write(comment)

def update_mesh():
    """
    """
    emitter_lines = open("mesh.txt").readlines()
    results_lines = open("results_data.00000001").readlines()
    results_lines = results_lines[results_lines.index("ASCII\n")+1:]
    remove_ids = [line.split()[2] for line in results_lines]
    new_emitter_lines = []
    i = 1
    for line in emitter_lines[1:-1]:
        sl = line.split()
        if i in remove_ids:
            sl = [sl[0], sl[1], sl[2], "0", sl[4]]
        new_line = "{}\n".format("	".join(sl))
        new_emitter_lines.append(new_line)
        i += 1
    with open("updated_mesh.txt", "w") as ue:
        ue.write(emitter_lines[0])
        for line in new_emitter_lines:
            ue.write(line)
        ue.write(emitter_lines[-1])


def write_meshgen_ini():
    with open("meshgen.ini", "w") as mgn:
        for line in MGN_INI_LINES:
            mgn.write(line)


def relax_emitter(n_nodes):

    print("Converting emitter to LAMMPS structure")
    convert_emitter_to_lammps("updated_mesh.txt", find_surface_atoms())
    print("Writing LAMMPS input file")
    write_lammps_input_file("data.emitter")

    print("Running LAMMPS")
    _ = subprocess.check_output([LAMMPS_CMD, "-l", "log.lammps",
                                 "-i", "in.emitter_relax"])
    print("Converting LAMMPS structure back to emitter")
    assign_ids_by_cn("relaxed_emitter.lmp")
    convert_lammps_to_emitter("relaxed_emitter.lmp", n_nodes)
    add_original_vacuum_nodes()
#    remove_duplicate_nodes("emitter_{}.txt".format(step_number))


def find_surface_atoms():
    surface_file = [f for f in os.listdir(os.getcwd()) if "surface_data" in f][-1]
    surface_lines = open(surface_file).readlines()
    surface_numbers = []
    for line in surface_lines[5:]:
        split_line = line.split()
        if len(split_line) == 0:
            break
        elif split_line[1] == "10":
            surface_numbers.append(split_line[0])
    return surface_numbers


def convert_emitter_to_lammps(emitter_file, surface_numbers):

    emitter_lines = open(emitter_file).readlines()
    atom_lines = [l for l in emitter_lines[1:-1] if l.split()[-2] not in
                  ["0", "1", "2", "3"]]
    atom_types = [t.split("=")[0][0] for t in emitter_lines[-1].split()[1:]]

    atom_coords = []
    n = 1
    for line in atom_lines:
        split_line = line.split()
        atom_coords.append(
            [str(n), split_line[-2][0]]+
            [str(float(x)*1e10) for x in split_line[:3]]
        )
        n += 1

    x_coords = [float(a[2]) for a in atom_coords]
    y_coords = [float(a[3]) for a in atom_coords]
    z_coords = [float(a[4]) for a in atom_coords]
    xlim = (min(x_coords)-10, max(x_coords)+10)
    ylim = (min(y_coords)-10, max(y_coords)+10)
    zlim = (min(z_coords)-10, max(z_coords)+10)

    with open("data.emitter", "w") as dat:
        dat.write("LAMMPS Emitter\n\n")
        dat.write("{} atoms\n\n".format(len(atom_lines)))
        dat.write("{} atom types\n\n".format(len(atom_types)))
        dat.write("{} {} xlo xhi\n".format(xlim[0], xlim[1]))
        dat.write("{} {} ylo yhi\n".format(ylim[0], ylim[1]))
        dat.write("{} {} zlo zhi\n\n".format(zlim[0], zlim[1]))
        dat.write("Masses\n\n")
        for elt in atom_types:
            dat.write("{} {}\n".format(elt, MASSES[elt]))
        dat.write("\nAtoms\n\n")
        fixed_indices = []
        for atom in atom_coords:
            if atom[0] not in surface_numbers:
                fixed_indices.append(atom[0])
            dat.write("{}\n".format(" ".join(atom)))
    with open("fixed_indices.txt", "w") as fi:
        for index in fixed_indices:
            fi.write("{}\n".format(index))


def convert_lammps_to_emitter(relaxed_structure_file, n_nodes):
    xyz_lines = open(relaxed_structure_file).readlines()
    with open("relaxed_emitter.txt", "w") as e:
        e.write("ASCII {} 0 0\n".format(n_nodes))
        i = 1
        for line in xyz_lines[9:]:
            sl = line.split()
            x = str(float(sl[0])*1e-10)
            y = str(float(sl[1])*1e-10)
            z = str(float(sl[2])*1e-10)
            atom_type = sl[3]
#            atom_id = sl[4]
            e.write("{}\n".format("	".join([x, y, z, atom_type])))
            i += 1


def add_original_vacuum_nodes():
    original_emitter_lines = open("../0/emitter.txt").readlines()
    original_vacuum_lines = [
        l for l in original_emitter_lines[1:] if l[0] == "#" or
        l.split()[3] in ["0", "2"]
    ]

    emitter_lines = open("relaxed_emitter.txt").readlines()
    emitter_lines = [
        l for l in emitter_lines[1:-1] if l.split()[3] not in
        ["0", "2"]
    ]
    n_nodes = len(emitter_lines)+len(original_vacuum_lines)-1

    with open("relaxed_emitter.txt", "w") as e:
        e.write("ASCII {} 0 0\n".format(n_nodes))
        for line in emitter_lines + original_vacuum_lines:
            e.write(line)


def remove_duplicate_nodes(emitter_file):
    e_lines = open(emitter_file).readlines()
    original_n_nodes = int(e_lines[0].split()[1])
    unique_coords, unique_lines = [], []
    n_duplicates = 0
    for line in e_lines:
        sl = line.split()
        if len(sl) > 2 and " ".join(sl[:3]) not in unique_coords:
            unique_coords.append(" ".join(sl[:3]))
            unique_lines.append(line)
        else:
            n_duplicates += 1
    with open(emitter_file, "w") as e:
        e.write("ASCII {} 0 0\n".format(original_n_nodes-n_duplicates))
        for line in unique_lines:
            e.write(line)
        e.write(e_lines[-1])

def write_lammps_input_file(structure_file):
    fixed_indices = [l.replace("\n", "") for l in open("fixed_indices.txt").readlines()]
    with open("in.emitter_relax", "w") as er:
        er.write("# Emitter Relaxation\n\n")
        er.write("units real\natom_style atomic\n\nread_data {}\n\n".format(structure_file))
        er.write("pair_style meam/c\n")
        er.write("pair_coeff * * /u/mashton/software/lammps/potentials/library.meam W NULL W\n\n")
        er.write("neighbor 1.0 bin\n")
        er.write("group inner id {}\n".format(" ".join([i for i in fixed_indices])))
        er.write("velocity inner set 0 0 0\n")
        er.write("fix frozen inner setforce 0 0 0\n\n")
        er.write("compute cnum all coord/atom cutoff 3.0\n")
        er.write("dump 1 all custom 1000 cnum.dump c_cnum\n")
        er.write("fix 1 all nve\n")
        er.write("minimize 1e-8 1e-8 1000 1000\n")
        er.write("compute 1 all coord/atom cutoff 4.0\n")
        er.write("write_dump all custom relaxed_emitter.lmp x y z type c_cnum")
