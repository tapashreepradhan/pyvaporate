import subprocess

import os

from pyvaporate.evaluate import assign_ids_by_cn
from pyvaporate.mgn import mgn_ini_lines

from monty.serialization import loadfn


def call_tapsim(node_file, setup):
    """
    Run the TAPSim program, starting with building a voronoi mesh for
    an emitter (node) file and then running a set number of
    evaporation steps before stopping.
    """

    write_meshgen_ini()

    _ = subprocess.check_output(
        [setup["evaporation"]["meshgen_bin"], node_file, "mesh.txt",
         "--create-config-template=mesh.cfg", "--write-ascii"]
    )

    with open("mesh.txt", "a") as m:
        comment = "#"
        for ID in setup["id_dict"]:
            comment += " {}={}".format(ID, setup["id_dict"][ID])
        m.write(comment)

    sm_lines = open("mesh.cfg").readlines()
    with open("mesh.cfg", "w") as sm:
        edit = False
        for line in sm_lines:
            if "ID" in line and line.split()[-1] not in ["0", "1", "2", "3"]:
                ID = line.split()[-1]
                BASE_ID = str(int(ID) - int(ID[-1]))
                elt = setup["id_dict"][BASE_ID]
                name = elt+"_"+ID[-1]
                edit = True
            if edit:
                if "***_SET_NAME_HERE_***" in line:
                    line = line.replace("***_SET_NAME_HERE_***", name)
                elif "***_SET_MASS_HERE_***" in line:
                    line = line.replace("***_SET_MASS_HERE_***", str(setup["emitter"]["elements"][elt]["mass"]))
                elif "EVAPORATION_CHARGE_STATE" in line:
                    line = line.replace("1", str(setup["emitter"]["elements"][elt]["charge"]))
                elif "***_SET_EVAPORATION_FIELD_STRENGTH_HERE_***" in line:
                    line = line.replace(
                        "***_SET_EVAPORATION_FIELD_STRENGTH_HERE_***",
                        setup["emitter"]["elements"][elt]["e_fields"][int(ID[-1])]
                    )
                sm.write(line)
            else:
                sm.write(line)

    _ = subprocess.check_output(
        [setup["evaporation"]["tapsim_bin"], "evaporation", "mesh.cfg", "mesh.txt",
         "--event-limit={}".format(setup["evaporation"]["events_per_step"]), "--write-ascii"]
    )
    update_mesh()


def update_mesh():
    """
    Remove evaporated nodes from the TAPSim mesh.
    """
    emitter_lines = open("mesh.txt").readlines()
    results_lines = open("results_data.00000001").readlines()
    results_lines = results_lines[results_lines.index("ASCII\n")+1:]
    remove_ids = [line.split()[2] for line in results_lines]
    new_emitter_lines = []
    i = 1
    for line in emitter_lines[1:-1]:
        sl = line.split()
        if str(i) in remove_ids:
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
    """
    Write the default meshgen.ini file so that
    meshgen doesn't ask for it mid-run and
    interrupt the python process.
    """
    with open("meshgen.ini", "w") as mgn:
        for line in mgn_ini_lines:
            mgn.write(line)


def call_lammps(n_nodes, setup):
    """
    Convert a TAPSim emitter node file to a LAMMPS
    structure (Only the actual atoms, not the vacumu nodes,
    etc) and then make a call to LAMMPS to relax it.
    """

    convert_emitter_to_lammps("updated_mesh.txt", find_surface_atoms(), setup)
    write_lammps_input_file("data.emitter", setup)

    _ = subprocess.check_output([setup["lammps"]["bin"], "-l", "log.lammps",
                                 "-i", "in.emitter_relax"])
    assign_ids_by_cn("relaxed_emitter.lmp")
    convert_lammps_to_emitter("relaxed_emitter.lmp", n_nodes)
    add_original_vacuum_nodes()


def find_surface_atoms():
    """
    Read from the TAPSim-generated surface_data file
    which atoms are at the emitter's surface. This is
    convenient if these are the only atoms you want to
    allow to relax in LAMMPS (surface_only: true in your
    setup.yaml) to speed up the relaxation.
    """

    surface_files = [f for f in os.listdir(os.getcwd()) if "surface_data" in f]
    if len(surface_files):
        surface_file = surface_files[-1]
        surface_lines = open(surface_file).readlines()
        surface_numbers = []
        for line in surface_lines[5:]:
            split_line = line.split()
            if len(split_line) == 0:
                break
            elif split_line[1] == "10":
                surface_numbers.append(split_line[0])
    else:
        surface_numbers = []
    return surface_numbers


def convert_emitter_to_lammps(emitter_file, surface_numbers, setup):
    """
    Convert a TAPSim emitter node file to a LAMMPS
    structure file.
    """

    emitter_lines = open(emitter_file).readlines()
    atom_lines = [l for l in emitter_lines[1:-1] if l.split()[-2] not in
                  ["0", "1", "2", "3"]]
    atom_types = [t.split("=")[0][0] for t in emitter_lines[-1].split()[1:]]
    atom_names = [t.split("=")[1] for t in emitter_lines[-1].split()[1:]]

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
        for i in range(len(atom_types)):
            dat.write("{} {}\n".format(
                atom_types[i],
                setup["emitter"]["elements"][atom_names[i]]["mass"])
            )
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
    """
    Convert a relaxed LAMMPS structure file back to
    a TAPSim emitter node file.
    """

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
    """
    Add the vaccuum, etc. nodes back to a TAPSim
    emitter node file created from a LAMMPS structure,
    which neither needs nor has these nodes.
    """
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


def write_lammps_input_file(structure_file, setup):
    """
    Write the input file specifying the type
    of relaxation to perform in LAMMPS.
    """

    etol = setup["lammps"]["minimize"]["etol"]
    ftol = setup["lammps"]["minimize"]["ftol"]
    maxiter = setup["lammps"]["minimize"]["maxiter"]
    maxeval = setup["lammps"]["minimize"]["maxeval"]
    pot = setup["lammps"]["potentials_location"]
    elts = " ".join(setup["emitter"]["elements"])

    fixed_indices = [l.replace("\n", "") for l in open("fixed_indices.txt").readlines()]
    with open("in.emitter_relax", "w") as er:
        er.write("# Emitter Relaxation\n\n")
        er.write("units real\natom_style atomic\n\nread_data {}\n\n".format(structure_file))
        er.write("pair_style meam/c\n")
        er.write("pair_coeff * * %s %s NULL %s\n\n" % (pot, elts, elts))
        er.write("neighbor 1.0 bin\n")
        if setup["lammps"]["minimize"]["surface_only"] == True:
            er.write("group inner id {}\n".format(" ".join([i for i in fixed_indices])))
            er.write("velocity inner set 0 0 0\n")
            er.write("fix frozen inner setforce 0 0 0\n\n")
        er.write("compute cnum all coord/atom cutoff 3.0\n")
        er.write("dump 1 all custom 1000 cnum.dump c_cnum\n")
        er.write("fix 1 all nve\n")
        er.write("minimize %s %s %s %s\n" % (etol, ftol, maxiter, maxeval))
        er.write("compute 1 all coord/atom cutoff 3.0\n")
        er.write("write_dump all custom relaxed_emitter.lmp x y z type c_cnum")
