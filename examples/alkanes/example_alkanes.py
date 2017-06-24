from matplotlib import pyplot as plt
from tau import tau

alkanes = ["methane.xyz", "ethane.xyz", "propane.xyz", "butane.xyz", "pentane.xyz"]
ccs = [tau.pa_ccs(xyzfile=alkane) for alkane in alkanes]

plt.plot([1, 2, 3, 4, 5], ccs, 'bo')
plt.xlabel("chain length"); plt.ylabel("CCS / A^2")
plt.savefig("example_alkanes.png")
plt.show()

