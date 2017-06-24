
This project aims to provide a simple to use python library for the calculation of
collision cross section (CCS) values for ion mobility spectrometry.
Currently, only the projection approximation (PA) method is implemented.
For example, the simple python script
```
from matplotlib import pyplot as plt
from tau import tau

alkanes = ["methane.xyz", "ethane.xyz", "propane.xyz", "butane.xyz", "pentane.xyz"]
ccs = [tau.pa_ccs(xyzfile=alkane) for alkane in alkanes]

plt.plot([1, 2, 3, 4, 5], ccs, 'bo')
plt.savefig("example_alkanes.png")
plt.show()
```
produces the output:


![alt text](https://raw.githubusercontent.com/Laangsam/tau/examples/alkanes/example_alkanes.png)


