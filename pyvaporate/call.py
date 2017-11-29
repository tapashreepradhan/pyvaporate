import subprocess

import os


LAMMPS_CMD = os.environ["LAMMPS_CMD"]
TAPSIM_CMD = os.environ["TAPSIM_BIN"] + "/tapsim"
MESHGEN_CMD = os.environ["TAPSIM_BIN"] + "/meshgen"
MGN_INI_LINES = open("{}/meshgen.ini".format(
    "/".join(__file__.split("/")[:-1]))
).readlines()


def call_tapsim(node_file, n_events):
    """
    Run the TAPSim program, starting with building a voronoi mesh for
    an emitter (node) file and then running a set number of
    evaporation steps before stopping.
    """

    write_meshgen_ini()


    _ = subprocess.check_output(
        [MESHGEN_CMD, node_file, "sampleMesh.bin",
         "--create-config-template=sampleMesh.cfg"]
    )

    sm_lines = open("sampleMesh.cfg").readlines()
    with open("sampleMesh.cfg", "w") as sm:
        #TODO: Actually edit the cfg file with atom names and information.
        for line in sm_lines:
            sm.write(line)

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
