<img src="https://s3.amazonaws.com/mashton/pyvaporate_logo.svg" style="position: relative; width: 50%; left: 25%"/>

Python Package for simulating full-scale (DFT-MD-MC) field evaporation.

#Installation
------
```bash
$ git clone https://github.com/ashtonmv/pyvaporate
```

Add pyvaporate to your `$PYTHONPATH`.

Pyvaporate requires the installation of LAMMPS (molecular dynamics) and TAPSim
(monte carlo field evaporation).

##LAMMPS Installation##

```
$ git clone https://github.com/lammps/lammps.git
$ cd lammps/src/
$ make yes-user-meamc
$ make mpi
```

Export the location of the create LAMMPS executable (`lmp_mpi`) to an
environment variable called `LAMMPS_CMD`.

##TAPSim Installation##

Visit (the University of Stuttgart's webpage)[http://www.uni-stuttgart.de/imw/mp/forschung/atom_probe_RD_center/software.en.html]
and download the TAPSim tarball. Unpack it and then execute the following:

```
$ cd tapsim_v1.0b_r3225/
$ make
```

Export the location of the `tapsim` and `meshgen` binaries created in `./bin`
to environment variables called `TAPSIM_CMD` and `MESHGEN_CMD`, respectively.

#Usage
------

Create a `.yaml` file with specifications of the system you wish to evaporate.
Default values exist for most available fields, but the user should familiarize
his/herself with these fields anyway. Below is a sample input file with all
the default values listed:

```
elements:
  W:
    percent_occ: 100  # Element symbol and percent of all sites occupied by that element
    e_fields:
      all: auto  # If "auto" is specified, all evaporation fields are looked up
                 # from PyVaporate's tables.
      # Otherwise, specify the evaporation field (N/m) for 
      # atoms of each coordination.
      1: auto
      2: auto
      3: auto
      4: auto
      5: auto
      6: auto
emitter:
  orientation:
    z: (1, 1, 0)
    y: auto
    x: auto
  radius: 100  # In Angstroms
  side_height: 50
evaporation:
  total_events: 100%  # total_events can be a percentage or an absolute number
  events_per_step: 10%  # Same goes for events_per_step
lammps:
  read_file: none  # Specify the path to a LAMMPS input file to use as a template for all MD relaxations. If not "none", this overrides the other commands in this section.
  potentials_location: /path/to/your library.meam  # A file in the lammps/potentials directory.
  minimize:
    etol: 1e-8
    ftol: 1e-8
    maxiter: 1000
    maxeval: 1000
```

After the input file is created, PyVaporate can be called simply as

```bash
$ pyvaporate input.yaml
```
