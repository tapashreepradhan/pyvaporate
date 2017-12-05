


def assign_ids_by_cn(lammps_file, id_dict):
    """
    Assign ID's to distinguish between atoms of various coordinations.
    """

    lammps_lines = open(lammps_file).readlines()
    with open(lammps_file, "w") as l:
        for line in lammps_lines[:9]:
            l.write(line)
        for line in lammps_lines[9:]:
            sl = line.split()
            original_atom_base = sl[-2]+"0"
            cn = int(sl[-1])
            if cn > 9:
                cn = 9
            atom_type = str(int(original_atom_base)+cn)
            sl[3] = atom_type
            l.write(" ".join(sl))
            l.write("\n")
