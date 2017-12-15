SETUP = {
    "basis": "BCC",
    "elements": {
        "W": {
            "percent_occ": 100,
            "charge": 3,
            "mass": 183.85,
            "e_fields": {
                0: 57e-9,
                1: 27e-9,
                2: 37e-9,
                3: 47e-9,
                4: 57e-9,
                5: 67e-9,
                6: 77e-9,
                7: 87e-9,
                8: 97e-9,
                9: 107e-9
            }
        }
    },
    "emitter": {
        "orientation": {
            "z": [1,1,0], "y": "auto", "x": "auto"
        },
        "radius": 100,
        "side_height": 50
    },
    "evaporation": {
        "tapsim_bin": "~/bin/tapsim",
        "meshgen_bin": "~/bin/meshgen",
        "total_events": "100%",
        "events_per_step": "10%"
    },
    "lammps": {
        "bin": "~/bin/lmp",
        "potentials_location": "~/software/lammps/potentials/library.meam",
        "input_template": "none",
        "minimize": {
            "etol": 1e-8, "ftol": 1e-8, "maxiter": 1000, "maxeval": 1000
        }
    }
}
