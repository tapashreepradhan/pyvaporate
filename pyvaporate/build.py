from ase.lattice.cubic import SimpleCubic, FaceCenteredCubic, BodyCenteredCubic

import sys
import math
import numpy as np
from random import randint


ELTS = {"10": "W"}
MASSES = {"W": "183.85"}


def build_emitter(element, basis, z_axis, filename="emitter.txt", x_axis="auto",
                  y_axis="auto", emitter_radius=100, emitter_side_height=50,
                  vacuum_radius=25, alloy=False):
    """
    Build an emitter (set of nodes, TAPSim style) based on an
    element, basis, orientation, and dimensions.
    """

    IDS = ["10"]

    if x_axis == "auto" and y_axis == "auto":
        if 0.99 < np.dot(z_axis, (1, 0, 0)) < 1.01:
            x_axis = tuple(np.cross(z_axis, (1, 0, 0)))
        else:
            x_axis = tuple(np.cross(z_axis, (0, 1, 0)))
        y_axis = tuple(np.cross(z_axis, x_axis))
        print(x_axis, y_axis, z_axis)

    R = emitter_radius + vacuum_radius

    if basis.lower() == "bcc":
        lattice = BodyCenteredCubic(
            size=(1, 1, 1),
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )
        pts = lattice.get_positions()
        supercell_dimensions = (
            int(math.ceil((2*R)/max([pt[0] for pt in pts])))+1,
            int(math.ceil((2*R)/max([pt[1] for pt in pts])))+1,
            int(math.ceil((emitter_side_height+R)/max([pt[2] for pt in pts])))+1
        )
        lattice = BodyCenteredCubic(
            size=supercell_dimensions,
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )
    elif basis.lower() == "fcc":
        lattice = FaceCenteredCubic(
            size=(1, 1, 1),
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )
        pts = lattice.get_positions()
        supercell_dimensions = (
            int(math.ceil((2*R)/max([pt[0] for pt in pts])))+1,
            int(math.ceil((2*R)/max([pt[1] for pt in pts])))+1,
            int(math.ceil((emitter_side_height+R)/max([pt[2] for pt in pts])))+1
        )
        lattice = FaceCenteredCubic(
            size=supercell_dimensions,
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )
    elif basis.lower() == "sc":
        lattice = SimpleCubic(
            size=(1, 1, 1),
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )
        pts = lattice.get_positions()
        supercell_dimensions = (
            int(math.ceil((2*R)/max([pt[0] for pt in pts])))+1,
            int(math.ceil((2*R)/max([pt[1] for pt in pts])))+1,
            int(math.ceil((emitter_side_height+R)/max([pt[2] for pt in pts])))+1
        )
        lattice = SimpleCubic(
            size=supercell_dimensions,
            directions=[x_axis, y_axis, z_axis],
            symbol=element
        )

    pts = lattice.get_positions()

    cx = np.mean([pt[0] for pt in pts])
    cy = np.mean([pt[1] for pt in pts])
    min_z = min([pt[2] for pt in pts])


    emitter_points, vacuum_points, bottom_points = [], [], []
    number = 0
    for pt in pts:
        pt = [pt[0], pt[1], pt[2]-min_z]
        if (pt[2] < 1e-5 and (pt[0]-cx)**2+(pt[1]-cy)**2<R**2):
            pt = [pt[0], pt[1], 0.0]
            pt = [i*1e-10 for i in pt]
            pt.append(2)  # bottom ID
            pt.append(number)
            number += 1
            bottom_points.append(pt)
        elif (pt[2]<emitter_side_height and (pt[0]-cx)**2+(pt[1]-cy)**2\
                <emitter_radius**2):
            pt = [i*1e-10 for i in pt]
            pt.append(10)  # emitter ID
            pt.append(number)
            number += 1
            emitter_points.append(pt)
        elif (pt[0]-cx)**2+(pt[1]-cy)**2+(pt[2]-emitter_side_height)**2\
                <emitter_radius**2:
            pt = [i*1e-10 for i in pt]
            pt.append(10)  # emitter ID
            pt.append(number)
            number += 1
            emitter_points.append(pt)
        elif (pt[2]<emitter_side_height and (pt[0]-cx)**2+(pt[1]-cy)**2<R**2):
            pt = [i*1e-10 for i in pt]
            pt.append(0)  # vacuum ID
            pt.append(number)
            number += 1
            vacuum_points.append(pt)
        elif (pt[0]-cx)**2+(pt[1]-cy)**2+(pt[2]-emitter_side_height)**2<R**2:
            pt = [i*1e-10 for i in pt]
            pt.append(0)  # vacuum ID
            pt.append(number)
            number += 1
            vacuum_points.append(pt)

    if alloy:
        conc = 0.05  # 10%
        substitution_indices = []
        num_atoms = len(emitter_points)
        num_substitution_atoms = math.floor(conc*num_atoms)
        j = 0
        while j < num_substitution_atoms:
            x = randint(0, num_atoms)
            if x not in substitution_indices:
                substitution_indices.append(x)
                j += 1
        for i in substitution_indices:
            emitter_points[i][3] = 11

    with open(filename, "w") as e:
        n_nodes = number
        e.write("ASCII {} 1 0\n".format(n_nodes))
        for pt in emitter_points + vacuum_points + bottom_points:
                # It's required that the coordinates be
                # separated by a tab character (^I), not
                # by regular spaces.
                e.write("	".join([str(i) for i in pt]))
                e.write("\n")
        comment = ["#"]
        comment += ["{}={}".format(ID, ELTS[ID]) for ID in IDS]
e.write("{}\n".format(" ".join(comment)))
