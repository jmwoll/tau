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

import unittest
import math
from tau import tau


class TauTest(unittest.TestCase):


    def test_pa_ccs_atom(self):
        ccs_val = tau.pa_ccs(xyzstring="C 0 0 0",radii='mobcal')
        expected = math.pi*(tau.parameters['mobcal']['C'])**2
        print(ccs_val, 'should be', expected)
        self.assertTrue(abs(ccs_val - expected) < 0.9)
        ccs_val = tau.pa_ccs(xyzstring="C 1 1 1",radii='mobcal')
        print(ccs_val)
        self.assertTrue(abs(ccs_val - expected) < 0.9)





if __name__ == '__main__':
    unittest.main()



