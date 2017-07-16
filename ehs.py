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

# Program for calculating the CCS via the exact hard sphere scattering model.
from tau import tau
import random
import math

class Line(object):

    def __init__(self, origin=None, direction=None):
        self.origin=origin; self.direction=direction

    def __str__(self):
        return ('<Line origin: {} | direction: {}>'.format(self.origin, self.direction))

class Sphere(object):

    def __init__(self, center=None, radius=None):
        self.radius=radius; self.center=center

# EHS parameters for the elements
parameters = {
'mobcal': {
    #'H': 2.2, 'C': 2.7, 'N': 2.7, 'O': 2.7 # PA values 120.9
    'H': 2.38, 'C': 3.043, 'N': 3.043, 'O': 3.043
}
}


def ehs_ccs(xyzfile=None, xyzstring=None, radii=None,
    num_rotamers=500):
    global parameters
    print('WARNING: EHS has not been correctly implemented (yet)!')
    assert xyzfile is not None or xyzstring is not None, (
        "the xyzfile must be specified")
    if radii is None or radii.lower() == 'siu_guo_2010':
        tau.log("No radii specified loading default radii:")
        tau.log("cite these parameters as:")
        #tau.log("Chi-Kit Siu Yuzhu Guo, Irine S. Saminathan, Alan C. Hopkinson," +
        #    " and K. W. Michael Siu*J. Phys. Chem. B, Vol. 114, No. 2, 2010.")
        tau.log("parameters are:\n")
        #tau.log(parameters['siu_guo_2010'])
        radii = parameters['mobcal'] #['siu_guo_2010']
    else:
        radii = parameters[radii]

    mol = (tau.load_molecule(xyzfile=xyzfile) if xyzfile is not None else
        tau.load_molecule(xyzstring=xyzstring))
    ccs_sum = 0.0
    print(mol)
    for _ in range(num_rotamers):
        randrot = [4 * math.pi * random.random() for _ in range(3)]
        mol = tau.rotate(mol, rot_x=randrot[0], rot_y=randrot[1],
            rot_z=randrot[2])
        min_x = min([atm[1] for atm in mol]) - 5
        max_x = max([atm[1] for atm in mol]) + 5
        min_y = min([atm[2] for atm in mol]) - 5
        max_y = max([atm[2] for atm in mol]) + 5
        ccs_sum += ehs_ccs_rotamer(mol, min_x, max_x, min_y, max_y, radii=radii)
    return ccs_sum / float(num_rotamers)


def fast_ehs_ccs_rotamer(molecule, min_x, max_x, min_y, max_y, trials=5000, radii=None):
    rand_xs = np.random.uniform(min_x, max_x, size=trials)
    rand_ys = np.random.uniform(min_y, max_y, size=trials)
    min_z = min([atom[3] for atom in molecule]) - 20
    n = [0, 0, 1]

    ray_dir = np.full((trials, 3), n)
    # ray_org = rand_x, rand_y, min_z
    # now find the first intersection of this ray
    # with all atoms in the molecule
    


def ehs_ccs_rotamer(molecule, min_x, max_x, min_y, max_y, trials=5000, radii=None):
    assert(radii is not None)
    max_min_x = max_x - min_x
    max_min_y = max_y - min_y
    min_z = min([atom[3] for atom in molecule]) - 20
    n = (0, 0, 1)
    ccs_sum = 0.0
    for _ in range(trials):
        rand_x = random.random() * max_min_x + min_x
        rand_y = random.random() * max_min_y + min_y
        ccs_sum += 2 * (1 - math.cos(scatter_angle(Line(direction=n,origin=(rand_x,rand_y,min_z)), molecule, radii=radii)))
    return (ccs_sum / float(trials)) * (max_x-min_x) * (max_y-min_y)


def scatter_angle(ray, molecule, radii=None):
    assert(radii is not None)
    scattered_ray = last_ray(ray, molecule, radii=radii)
    a_x, a_y, a_z = ray.direction; b_x, b_y, b_z = scattered_ray.direction
    ab = a_x * b_x + a_y * b_y + a_z * b_z
    abs_a = (a_x * a_x + a_y * a_y + a_z * a_z)**0.5
    abs_b = (b_x * b_x + b_y * b_y + b_z * b_z)**0.5
    return math.acos(ab / (abs_a * abs_b))


def last_ray(ray, molecule, radii=None, max_iterations=8):
    assert(radii is not None)
    iter_count = 0
    while True:
        new_ray_dir, new_ray_origin = ray_trace(ray, molecule, radii=radii)

        if new_ray_dir is None:
            return ray

        #new_ray_origin = (new_ray_scalar * ray.direction[0] + ray.origin[0],
        #    new_ray_scalar * ray.direction[1] + ray.origin[1],
        #    new_ray_scalar * ray.direction[2] + ray.origin[2])
        new_ray_len = (new_ray_dir[0] ** 2 + new_ray_dir[1] ** 2 + new_ray_dir[2]**2)**0.5
        new_ray_dir = (new_ray_dir[0] / new_ray_len, new_ray_dir[1] / new_ray_len,
            new_ray_dir[2] / new_ray_len)
        ray = Line(direction=new_ray_dir, origin=new_ray_origin)

        iter_count += 1
        if iter_count >= max_iterations:
            return ray




def ray_trace(ray, molecule, radii=None):
    assert(radii is not None)
    coll_scalar, coll_atom = first_collision(ray, molecule, radii=radii)
    if coll_scalar is None:
        return None, None
    # " coll_point = ray.direction * d + ray.origin "
    coll_point = (ray.direction[0]*coll_scalar + ray.origin[0],
                  ray.direction[1]*coll_scalar + ray.origin[1],
                  ray.direction[2]*coll_scalar + ray.origin[2])
    # the following code snippet calculates the reflected ray, w,
    # given the input line v,
    # (see https://math.stackexchange.com/questions/2334939/reflection-of-line-on-a-sphere/2334963?noredirect=1#comment4807112_2334963)
    # w = 2 * [(v * (x-c))] * (x-c) -v
    #         ^           ^ scalar product
    # lets first evaluate the scalar part
    c_to_x = (coll_point[0]-coll_atom[1],coll_point[1]-coll_atom[2],
           coll_point[2]-coll_atom[3])
    v = ray.direction
    scalar_fac = 2 * (v[0] * c_to_x[0] + v[1] * c_to_x[1] + v[2] * c_to_x[2])
    return (scalar_fac*c_to_x[0] - v[0], scalar_fac*c_to_x[1] - v[1],
            scalar_fac*c_to_x[2] - v[2]), coll_point


def first_collision(ray, molecule, radii=None):
    assert(radii is not None)
    # finds the first collision of the ray with the molecule
    atm_coll_dct = {}
    for atm in molecule:
        intscts = line_sphere_intersections(ray, Sphere(center=(atm[1],atm[2],atm[3]), radius=radii[atm[0]]))
        # Only intersection in the ray direction, not backwards
        if intscts:
            intscts = [intsct for intsct in intscts if intsct > 0]
        if intscts:
            atm_coll_dct[min(intscts)] = atm
    if atm_coll_dct:
        # only accept collisions that don't turn backwards,
        # as the buffer gas only moves in one direction.
        # otherwise, it might bounce backwards
        # artificially between atoms.
        #positive_keys = [key for key in atm_coll_dct.keys() if key > 0]
        #if positive_keys:
        #fst_coll = min(positive_keys)
        assert(len([key for key in atm_coll_dct.keys() if key <= 0]) == 0)
        fst_coll = min(atm_coll_dct.keys())
        return fst_coll, atm_coll_dct[fst_coll]
    return None, None


def line_sphere_intersections(line, sphere):
    #print("------------->\n",line)
    #print('<---------->')
    oc = (line.origin[0]-sphere.center[0],
        line.origin[1]-sphere.center[1], line.origin[2]-sphere.center[2])
    oc_sum = oc[0] + oc[1] + oc[2]
    oc_abs_sq = (oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2])
    oc_sq = oc_sum**2
    radicant = oc_sq - oc_abs_sq + sphere.radius**2 # is this alright?
    if radicant < 0:
        return False
    return -oc_sum + (radicant)**0.5, -oc_sum - (radicant)**0.5
