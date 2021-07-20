import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_csv("diag.csv",header=0)

e = np.sqrt((d["rms_u"]**2 + d["rms_v"]**2 + d["rms_w"]**2) / 4.0 + 10.0 * d["rms_cmt"]**2 / 4.0 / 1.0)

start = np.argmax(d["time"] > 20.0)
end = np.argmax(d["time"] > 24.0)
period_start = np.argmax(e[start:end])

start = np.argmax(d["time"] > 32.0)
end = np.argmax(d["time"] > 36.0)
period_end = np.argmax(e[start:end])

period = d["time"][period_end] - d["time"][period_start]
max_cycles = np.floor((100.0 - d["time"][period_start]) / period)
start = period_start
end = np.argmax(d["time"] > (d["time"][period_start] + period * max_cycles))

print(start,end,period)

p = np.polyfit(d["time"][start:end],np.log(e[start:end]),1)
print(p)

plt.scatter([d["time"][period_start],d["time"][period_end]],[e[period_start],e[period_end]])

plt.plot(d["time"],e)
plt.plot(d["time"],np.exp(np.polyval(p,d["time"])))
plt.yscale("log")
plt.show()
