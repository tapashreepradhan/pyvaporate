from monty.serialization import loadfn

from pyvaporate.build import build_emitter
from pyvaporate.call import call_tapsim, call_lammps
from pyvaporate import SETUP

import os

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
    if len(elements) > 1:
        for e in elements[1:]:
            alloy[e] = SETUP["emitter"]["elements"][e]["fract_occ"]

    step_number = 0
    if "%" in n_events_total:
        total_percent = float(n_events_total.replace("%",""))/100.
        n_events_total = np.inf
    if "%" in n_events_per_step:
        step_percent = float(n_events_per_step.replace("%",""))/100.
        n_events_per_step = 0
    while step_number * SETUP["evaporation"]["events_per_step"] < SETUP["evaporation"]["total_events"]:
        if not os.path.isdir(str(step_number)):
            os.mkdir(str(step_number))
        os.chdir(str(step_number))
        print("\nSTEP {}\n------".format(step_number))
        if step_number == 0:
            if SETUP["emitter"]["file"] != "none":
                print("Importing emitter from {}".format(SETUP["emitter"]["file"]))
                os.system("cp {} ./emitter.txt".format(SETUP["emitter"]["file"]))
            else:
                print("Building initial emitter")
                basis = SETUP["emitter"]["basis"]
                emitter_radius = SETUP["emitter"]["radius"]
                emitter_side_height = SETUP["emitter"]["side_height"]
                z_axis = SETUP["emitter"]["orientation"]["z"]
                y_axis = SETUP["emitter"]["orientation"]["y"]
                x_axis = SETUP["emitter"]["orientation"]["x"]
                build_emitter(
                    element=elements[0], basis=basis, z_axis=z_axis,
                    filename="emitter.txt", emitter_radius=emitter_radius,
                    emitter_side_height=emitter_side_height, alloy=alloy
                )
            n_atoms = len([l for l in open("emitter.txt").readlines()[1:-1] if
                           l.split()[3] not in ["0", "1", "2", "3"]])
            if n_events_total == np.inf:
                SETUP["evaporation"]["total_events"] = math.ceil(total_percent * n_atoms)
            if n_events_per_step == 0:
                SETUP["evaporation"]["events_per_step"] = math.ceil(step_percent * n_atoms)
        else:
            os.system("cp ../{}/relaxed_emitter.txt emitter.txt".format(step_number-1))
        call_tapsim("emitter.txt", SETUP)
        n_atoms = int(open("emitter.txt").readlines()[0].split()[1])-n_events_per_step
        call_lammps(n_atoms, SETUP)
        step_number += 1
        os.chdir("../")
