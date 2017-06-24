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
from tau import tau

alkanes = ["methane.xyz", "ethane.xyz", "propane.xyz", "butane.xyz", "pentane.xyz"]
ccs = [tau.pa_ccs(xyzfile=alkane) for alkane in alkanes]

plt.plot([1, 2, 3, 4, 5], ccs, 'bo')
plt.xlabel("chain length"); plt.ylabel("CCS / A^2")
plt.savefig("example_alkanes.png")
plt.show()

