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
from timeit import default_timer as timer
import unittest
import math
from tau import tau



class TestPerformanceAlkanes(unittest.TestCase):


    def test_performance_pa_ccs_alkanes(self):
        alkanes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200]
        #reference_ccs = [27.499, 35.806, 42.457, 50.114, 57.079]
        ccs = []
        times = []
        for alkane in alkanes:
            start = timer()
            ccs.append(tau.pa_ccs(xyzfile="alkane_{}.xyz".format(alkane)))
            end = timer()
            times.append((end - start))


        #print(reference_ccs)
        plt.plot(alkanes, times, 'bo')
        plt.xlabel('chain length')
        plt.ylabel('time / s')
        plt.savefig('pa_ccs_time.svg')
        plt.show()

if __name__ == '__main__':
    unittest.main()