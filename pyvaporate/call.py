import subprocess

import os


LAMMPS_CMD = os.environ["LAMMPS_CMD"]
TAPSIM_CMD = os.environ["TAPSIM_BIN"] + "/tapsim"
MESHGEN_CMD = os.environ["TAPSIM_BIN"] + "/meshgen"
MGN_INI_LINES = open("{}/meshgen.ini".format(__file__.split("/")[:-1])).readlines()


def call_tapsim(node_file, n_events):
    """
    Run the TAPSim program, starting with building a voronoi mesh for
    an emitter (node) file and then running a set number of
    evaporation steps before stopping.
    """

    write_meshgen_ini()

    _ = subprocess.check_output(
        "{} {} sampleMesh.bin --create-config-template=sampleMesh.cfg".format(
            MESHGEN_CMD, node_file
        )
    )

    sm_lines = open("sampleMesh.cfg").readlines()
    with open("sampleMesh.cfg", "w") as sm:
        #TODO: Actually edit the cfg file with atom names and information.
        for line in sm_lines:
            sm.write(line)

    _ = subprocess.check_output(
        "{} evaporation sampleMesh.cfg sampleMesh.bin --event-limit={}".format(
            TAPSIM_CMD, n_events
        )
    )


def write_meshgen_ini():
    with open("meshgen.ini", "w") as mgn:
        for line in MGN_INI_LINES:
            mgn.write(line)


def call_lammps(emitter_file):

    #TODO: Build emitter_relax.in

    _ = subprocess.check_output("{} emitter_relax.in".format(LAMMPS_CMD))
