SETUP = {
    "basis": "BCC",
    "elements": {
        "W": {
            "percent_occ": 100,
            "e_fields": {
                "all": "auto"
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
        "potentials_location": "~/software/lammps/potentials",
        "input_template": "none",
        "minimize": {
            "etol": 1e-8, "ftol": 1e-8, "maxiter": 1000, "maxeval": 1000
        }
    }
}
