
This project aims to provide a simple to use python library for the calculation of
collision cross section (CCS) values for ion mobility spectrometry.
Currently, only the projection approximation (PA) method is implemented.

Please consider using https://github.com/jmwoll/goccs instead, which is a more up-to-date version of the project.

For example, the simple python script
```
from matplotlib import pyplot as plt
from tau import tau

alkanes = ["methane.xyz", "ethane.xyz", "propane.xyz", "butane.xyz", "pentane.xyz"]
ccs = [tau.pa_ccs(xyzfile=alkane) for alkane in alkanes]

plt.plot([1, 2, 3, 4, 5], ccs, 'bo')
plt.xlabel("chain length"); plt.ylabel("CCS / A^2")
plt.savefig("example_alkanes.png")
plt.show()
```
produces the output:
![Plot of alkane CCS values](https://github.com/Laangsam/tau/blob/master/examples/alkanes/example_alkanes.png)

where the CCS values are plotted on the vertical axis, and the chain length is the horizontal axis.
It was assumed here, that the xyz files are present in the script directory.

