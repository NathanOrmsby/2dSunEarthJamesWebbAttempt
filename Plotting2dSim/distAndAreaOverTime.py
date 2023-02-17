from matplotlib import pyplot as plt
import pandas as pd

# Desktop
f = "D:\\2dSunEarthChaos\\singleProbeData\\distArea.csv"

# Laptop
# fname = "C:\\Users\\natha\\eclipse-workspace\\3dSim\\unitCircle.csv"

df = pd.read_csv(f)
# print(df.head())

# Distance
y0 = df.loc[:, "dist"]

# Area
y1 = df.loc[:, "area"]

x = [i for i in range(len(y0))]

# Sphere Plot

# Distance Plot

plt.plot(x[:28184], y0[:28184])
plt.show()


# Area Plot

plt.plot(x, y0)
plt.show()