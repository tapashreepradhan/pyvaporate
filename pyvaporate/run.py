from monty.serialization import loadfn

from pyvaporate.build import build_emitter_from_scratch, build_emitter_from_file
from pyvaporate.call import call_tapsim, call_lammps
from pyvaporate import SETUP

import os
import sys
sys.stdout = open("pyvaporate.log", "w")

import numpy as np
import math


def yaml_run(config_file):
    """
    """

    CONFIGURATION = loadfn(config_file)
    for key in SETUP:
        if key in CONFIGURATION:
            SETUP[key] = CONFIGURATION[key]

    lammps = SETUP["lammps"]["bin"]

    tapsim = SETUP["evaporation"]["meshgen_bin"]
    meshgen = SETUP["evaporation"]["tapsim_bin"]
    n_events_total = SETUP["evaporation"]["total_events"]
    n_events_per_step = SETUP["evaporation"]["events_per_step"]

    elements = [e for e in SETUP["emitter"]["elements"]]
    alloy = {}
    SETUP["id_dict"] = {}
    n = 10
    for e in elements:
        SETUP["id_dict"][str(n)] = e
        n += 10
    if len(elements) > 1:
        for e in elements[1:]:
            alloy[e] = SETUP["emitter"]["elements"][e]["fract_occ"]

    step_number = 0
    if "%" in str(n_events_total):
        total_percent = float(n_events_total.replace("%",""))/100.
        n_events_total = np.inf
    if "%" in str(n_events_per_step):
        step_percent = float(n_events_per_step.replace("%",""))/100.
        n_events_per_step = 0
    while step_number * SETUP["evaporation"]["events_per_step"] <= SETUP["evaporation"]["total_events"]:
        if not os.path.isdir(str(step_number)):
            os.mkdir(str(step_number))
        os.chdir(str(step_number))
        print("\nSTEP {}\n------".format(step_number))
        if step_number == 0:
            if SETUP["emitter"]["source"]["node_file"] == "none" and SETUP["emitter"]["source"]["uc_file"] == "none":
                print("Building initial emitter")
                basis = SETUP["emitter"]["basis"]
                emitter_radius = SETUP["emitter"]["radius"]
                emitter_side_height = SETUP["emitter"]["side_height"]
                z_axis = SETUP["emitter"]["orientation"]["z"]
                y_axis = SETUP["emitter"]["orientation"]["y"]
                x_axis = SETUP["emitter"]["orientation"]["x"]
                build_emitter_from_scratch(
                    element=elements[0], basis=basis, z_axis=z_axis,
                    filename="emitter.txt", emitter_radius=emitter_radius,
                    emitter_side_height=emitter_side_height, alloy=alloy
                )

            elif SETUP["emitter"]["source"]["node_file"] != "none":
                print("Importing emitter from {}".format(SETUP["emitter"]["source"]["node_file"]))
                os.system("cp {} ./emitter.txt".format(SETUP["emitter"]["source"]["node_file"]))

            elif SETUP["emitter"]["source"]["uc_file"] != "none":
                print("Building emitter based on {}".format(SETUP["emitter"]["source"]["uc_file"]))
                emitter_radius = SETUP["emitter"]["radius"]
                emitter_side_height = SETUP["emitter"]["side_height"]
                z_axis = SETUP["emitter"]["orientation"]["z"]
                y_axis = SETUP["emitter"]["orientation"]["y"]
                x_axis = SETUP["emitter"]["orientation"]["x"]
                build_emitter_from_file(
                    SETUP["emitter"]["source"]["uc_file"], z_axis=z_axis,
                    filename="emitter.txt", emitter_radius=emitter_radius,
                    emitter_side_height=emitter_side_height
                )
            n_atoms = len([l for l in open("emitter.txt").readlines()[1:-1] if
                           l.split()[3] not in ["0", "1", "2", "3"]])
            if n_events_total == np.inf:
                SETUP["evaporation"]["total_events"] = math.ceil(total_percent * n_atoms)
            if n_events_per_step == 0:
                SETUP["evaporation"]["events_per_step"] = math.ceil(step_percent * n_atoms)
            lines = open("emitter.txt").readlines()
            with open("updated_mesh.txt", "w") as f:
                f.write(lines[0])
                for l in lines[1:-1]:
                    split_line = l.split()
                    split_line.append("0\n")
                    f.write(" ".join(split_line))
                f.write(lines[-1])
            call_lammps(n_atoms, SETUP)
        else:
            os.system("cp ../{}/relaxed_emitter.txt emitter.txt".format(step_number-1))
            call_tapsim("emitter.txt", SETUP)
            n_atoms = int(open("emitter.txt").readlines()[0].split()[1])-n_events_per_step
            call_lammps(n_atoms, SETUP)
            if SETUP["cleanup"] == True:
                os.system("rm trajectory_data.*")
                os.system("rm dump.*")
                os.system("rm dump")
                os.system("rm geometry.dat")

        step_number += 1
        os.chdir("../")
