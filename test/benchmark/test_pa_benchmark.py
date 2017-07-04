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

from matplotlib import pyplot as plt
import unittest
import math
import numpy as np
from tau import tau


def rms_error(predictions, targets):
    return np.sqrt(((np.array(predictions) - np.array(targets)) ** 2).mean())


class TestAlkanes(unittest.TestCase):


    def test_pa_ccs_alkanes(self):
        pa_ccs_ref = {
            "methane": 27.499, "ethane": 35.806, "propane": 42.457,
            "butane": 50.114, "pentane": 57.079, "adamantane": 64.799,
            "penguinone": 73.598, "barrelene": 56.733,
            "hirsutene": 88.677, "testosterone": 108.71,
            "paclitaxel": 231.00, "abacavir": 117.27,
            "ciclosporin": 286.18, "codeine": 104.30,
            "talinolol": 165.93, "yangonin": 115.52,
            "vernakalant": 145.82, "octabenzone": 149.15,
            "melatonin": 108.92, "epirubicin": 172.70,
            "alprenolol": 115.17, "trioxifene": 171.34,
            "methylaminoethanol": 53.907,
            "nanokid": 225.19,

        }

        pa_ccs_tau = {}
        for name in pa_ccs_ref:
            pa_ccs_tau[name] = tau.pa_ccs(xyzfile="{}.xyz".format(name),radii='mobcal')
        print(pa_ccs_tau)
        print(pa_ccs_ref)

        plt.plot([pa_ccs_tau[key] for key in sorted(pa_ccs_tau)], [pa_ccs_ref[key] for key in sorted(pa_ccs_ref)], 'bo')
        for key in pa_ccs_tau:
            px, py = pa_ccs_tau[key], pa_ccs_ref[key]
            plt.text(px, py+1, key, horizontalalignment='center')
        plt.savefig('benchmark.png')
        plt.show()

        error = rms_error([pa_ccs_tau[key] for key in sorted(pa_ccs_tau)], [pa_ccs_ref[key] for key in sorted(pa_ccs_ref)])
        mean_ccs = sum([pa_ccs_tau[key] for key in sorted(pa_ccs_tau)])/float(len(pa_ccs_tau))
        print('RMS error: {} ( {} % )'.format(error, 100*error/mean_ccs))
        self.assertTrue(error < 10, 'error must be below 10 angstrom^2')
        self.assertTrue(100*error/mean_ccs < 2, 'relative error must be below 2 %')


if __name__ == '__main__':
    unittest.main()
