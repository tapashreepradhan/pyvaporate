from monty.serialization import loadfn

from pyvaporate.build import build_emitter
from pyvaporate.call import call_tapsim, call_lammps
from pyvaporate import SETUP


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
    emitter_radius = SETUP["emitter"]["radius"]
    emitter_side_height = SETUP["emitter"]["side_height"]
    elements = [e for e in SETUP["elements"]]
    basis = SETUP["basis"]

    step_number = 0
    while step_number * n_events_per_step < n_events_total:
        if not os.path.isdir(str(step_number)):
            os.mkdir(str(step_number))
        os.chdir(str(step_number))
        print("\nSTEP {}\n------".format(step_number))
        if step_number == 0:
            print("Building initial emitter")
            build_emitter(
                element=elements[0], basis=basis, z_axis=z_axis,
                filename="emitter.txt", emitter_radius=emitter_radius,
                emitter_side_height=emitter_side_height
            )
        else:
            os.system("cp ../{}/relaxed_emitter.txt emitter.txt".format(step_number-1))
        call_tapsim("emitter.txt", n_events_per_step, id_dict={"10": elements[0]})
        n_atoms = int(open("emitter.txt").readlines()[0].split()[1])-n_events_per_step
        call_lammps(n_atoms)
        step_number += 1
        os.chdir("../")
