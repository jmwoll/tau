# Copyright (C) 2017  Jan Wollschläger <jmw.tau@gmail.com>
# This file is part of Tau.
#
# Tau is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import math
import random

_log_print = True


def log(msg):
    global _log_print
    if _log_print:
        print(msg)

numpy_loaded = True
try:
    import numpy as np
except:
    numpy_loaded = False
    log("Could not load numpy.")


parameters = {
# radii from mobcal
# J. Phys. Chem. 1996, 100, 16082-16086;
# J. Phys. Chem. A 1997, 101, 968.
# Chem. Phys. Lett. 1996, 261, 86-91.
'mobcal': {
    'H': 2.2, 'C': 2.7, 'N': 2.7, 'O': 2.7
},
# J. Phys. Chem. B, Vol. 114, No. 2, 2010.
'siu_guo_2010': {
    'H': 2.01, 'C': 2.35, 'N': 2.26, 'O': 2.26
}
}


def load_molecule(xyzfile=None, xyzstring=None):
    mol = []
    inp_str = None
    if xyzstring is None:
        with open(xyzfile, 'r') as fin:
            inp_str = fin.read()
    else:
        inp_str = xyzstring
    inp_str = inp_str.strip().replace('\t', ' ')
    while '  ' in inp_str:
        inp_str = inp_str.replace('  ', ' ')

    for i, lne in enumerate(inp_str.split('\n')):
        atm = lne.strip().split(' ')
        try:
            atm = (atm[0], float(atm[1]), float(atm[2]), float(atm[3]))
            mol.append(atm)
        except:
            if i not in [0,1]: log('skipping line:\n' + lne)
    return mol


def pa_ccs(xyzfile=None, xyzstring=None, radii=None,
    num_rotamers=5000, fast=None):
    global parameters
    if fast is None:
        fast = numpy_loaded
    assert xyzfile is not None or xyzstring is not None, (
        "the xyzfile must be specified")
    if radii is None or radii.lower() == 'siu_guo_2010':
        log("No radii specified loading default radii:")
        log("cite these parameters as:")
        log("Chi-Kit Siu Yuzhu Guo, Irine S. Saminathan, Alan C. Hopkinson," +
            " and K. W. Michael Siu*J. Phys. Chem. B, Vol. 114, No. 2, 2010.")
        log("parameters are:\n")
        log(parameters['siu_guo_2010'])
        radii = parameters['siu_guo_2010']
    else:
        radii = parameters[radii]

    mol = (load_molecule(xyzfile=xyzfile) if xyzfile is not None else
        load_molecule(xyzstring=xyzstring))
    ccs_sum = 0.0

    ccs_rotamer = fast_pa_ccs_rotamer if fast else pa_ccs_rotamer
    for _ in range(num_rotamers):
        randrot = [4 * math.pi * random.random() for _ in range(3)]
        mol = rotate(mol, rot_x=randrot[0], rot_y=randrot[1],
            rot_z=randrot[2])
        min_x = min([atm[1] for atm in mol]) - 5
        max_x = max([atm[1] for atm in mol]) + 5
        min_y = min([atm[2] for atm in mol]) - 5
        max_y = max([atm[2] for atm in mol]) + 5
        ccs_sum += ccs_rotamer(mol, radii, min_x, max_x, min_y, max_y)
    return ccs_sum / float(num_rotamers)


def rotate(mol, rot_x=None, rot_y=None, rot_z=None):
    rot_mol = []
    cx, cy, cz = 0, 0, 0
    for atm in mol:
        cx, cy, cz = cx + atm[1], cy + atm[2], cz + atm[3]

    cx, cy, cz = (cx / float(len(mol)),
        cy / float(len(mol)), cz / float(len(mol)))

    for atm in mol:
        s, x, y, z = atm[0], atm[1], atm[2], atm[3]
        x, y, z = x - cx, y - cy, z - cz
        y, z = (y * math.cos(rot_x) - z * math.sin(rot_x),
            y * math.sin(rot_x) + z * math.cos(rot_x))
        x, z = (x * math.cos(rot_y) + z * math.sin(rot_y),
            -x * math.sin(rot_y) + z * math.cos(rot_y))
        x, y = (x * math.cos(rot_z) - y * math.sin(rot_z),
            x * math.sin(rot_z) + y * math.cos(rot_z))
        x, y, z = x + cx, y + cy, z + cz
        rot_mol.append((s, x, y, z))
    assert(len(rot_mol) == len(mol))
    return rot_mol


def pa_ccs_rotamer(mol, radii, min_x, max_x, min_y, max_y, trials=5000):
    max_min_x = max_x - min_x
    max_min_y = max_y - min_y

    hits = 0
    hit_atom = False

    for _ in range(trials):
        rand_x = random.random() * max_min_x + min_x
        rand_y = random.random() * max_min_y + min_y
        hit_atom = False
        for atm in mol:
            #assert(not hit_atom)
            dx, dy = abs(rand_x - atm[1]), abs(rand_y - atm[2])
            if dx * dx + dy * dy < (radii[atm[0]]) ** 2:
                if not hit_atom:
                    hits += 1
                    hit_atom = True
                    break # no double/multiple hits

    return (hits / (float(trials))) * (max_x - min_x) * (max_y - min_y)


def fast_pa_ccs_rotamer(mol, radii, min_x, max_x, min_y, max_y, trials=10000):
     rand_xs = np.random.uniform(min_x, max_x, size=trials)
     rand_ys = np.random.uniform(min_y, max_y, size=trials)
     #deltas = (1 for randx,randy in zip(rand_xs, rand_ys) if any((abs(randx - atm[1])**2 + abs(randy - atm[2])**2 < (radii[atm[0]])**2 for atm in mol)))
     hits = np.zeros(trials)
     for atm in mol:
         atm_x, atm_y, atm_rad_sq = atm[1], atm[2], (radii[atm[0]] ** 2)
         # hits += len( ((rand_xs - atm_x)**2 + (rand_ys - atm_y)**2) < atm_rad_sq )
         hits = hits + ( ((rand_xs - atm_x)**2 + (rand_ys - atm_y)**2) < atm_rad_sq )

     return ( len([hit for hit in hits if hit > 0]) / float(trials) ) * (max_x - min_x) * (max_y - min_y)
