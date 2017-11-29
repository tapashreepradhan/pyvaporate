import subprocess

import os


LAMMPS_CMD = os.environ["LAMMPS_CMD"]
TAPSIM_CMD = os.environ["TAPSIM_BIN"] + "/tapsim"
MESHGEN_CMD = os.environ["TAPSIM_BIN"] + "/meshgen"
MGN_INI_LINES = open("{}/meshgen.ini".format(
    "/".join(__file__.split("/")[:-1]))
).readlines()
ELTS = {"10": "W"}
E_FIELDS = {"10": 57.1e-9}
MASSES = {"W": 183.85}
CHARGE_STATES = {"W": 3}

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
            else:
                sm.write(line)

    print("Running TAPSim")
    _ = subprocess.check_output(
        [TAPSIM_CMD, "evaporation", "sampleMesh.cfg", "sampleMesh.bin",
         "--event-limit={}".format(n_events)]
    )


def write_meshgen_ini():
    with open("meshgen.ini", "w") as mgn:
        for line in MGN_INI_LINES:
            mgn.write(line)


def call_lammps(emitter_file):

    #TODO: Build emitter_relax.in

    _ = subprocess.check_output([LAMMPS_CMD, "emitter_relax.in"])
