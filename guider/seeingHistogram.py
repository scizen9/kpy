import matplotlib.pyplot as plt
import numpy as np

filename = raw_input("File: ")
label = raw_input("Label: ")
f = open(filename)
fwhm = []
for line in f: 
    t = line.split('-')
    if t[0] != "nan":
        fwhm.append(float(t[0]))
plt.figure()
bins = np.arange(0.5, 3.0, 0.1)
plt.hist(fwhm, bins, label=label)
plt.xlabel("FWHM (arcseconds)")
plt.legend()
plt.show()
