import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FFMpegWriter, PillowWriter
import pandas as pd

# Desktop ffmpeg path
plt.rcParams['animation.ffmpeg_path'] = 'D:\\ffmpeg\\bin\\ffmpeg.exe'

metadata = dict(title='Movie', artist='me')
# writer = PillowWriter(fps=15, metadata=metadata)
writer = FFMpegWriter(fps=30, metadata=metadata)

fig = plt.figure()
l, = plt.plot([], [], 'k-')

# Read from csv files
# Desktop
f0 = "D:\\2dSunEarthChaos\\singleProbeData\\earthAndprobe.csv"
f1 = "D:\\2dSunEarthChaos\\singleProbeData\\perturbations.csv"
f2 = "D:\\2dSunEarthChaos\\singleProbeData\\gramSchmidtProbes.csv"
f3 = "D:\\2dSunEarthChaos\\singleProbeData\\distFromL2.csv"


df0 = pd.read_csv(f0)
df1 = pd.read_csv(f1)
df2 = pd.read_csv(f2)
df3 = pd.read_csv(f3)

# Earth
x0 = df0.loc[:, 'ex']
y0 = df0.loc[:, 'ey']

# Probe
# x1 = df0.loc[:, 'jwx']
# y1 = df0.loc[:, 'jwy']

# Perturbations
NUMPERTURBED = 5
x2s = [df1.loc[:, "x{}" .format(i)] for i in range(NUMPERTURBED)]
y2s = [df1.loc[:, "y{}" .format(i)] for i in range(NUMPERTURBED)]

# h2
h2v1x = df2.loc[:, "v1x"]
h2v1y = df2.loc[:, "v1y"]
h2v2x = df2.loc[:, "v2x"]
h2v2y = df2.loc[:, "v2y"]

# L2 point
x3 = df0.loc[:, 'x']
y3 = df0.loc[:, 'y']

# Point 1: Sun
x = 149597870700
y = 0

colors = ['g', 'r', 'c', 'm', 'y', 'k', 'black', 'aqua', 'lime', 'maroon', 'navy', 'olive', 'purple', 'grey', 'silver', 'teal', 'white', 'orange', 'hotpink', 'aquamarine']

# # # Point 1
# x1 = df1.loc[:, 'x']
# y1 = df1.loc[:, 'y']
# z1 = df1.loc[:, 'z']
# #
# # # Point 2
# x2 = df2.loc[:, 'x']
# y2 = df2.loc[:, 'y']
# z2 = df2.loc[:, 'z']

# Set limits
# ax.xlim(-10, 20)
# ax.ylim(-10, 20)
# ax.zlim(-10, 20)

with writer.saving(fig, "mp4s\\2dPaths.mp4", 200):
    for i in range(len(x0)):
        if i % 100 == 0:
            print(i)

        # Earth
        plt.plot(x0[:i + 1], y0[:i + 1], color="blue")
        plt.scatter(x0[i], y0[i], color="blue")

        # JW
        plt.plot(x1[:i + 1], y1[:i + 1], color="red")
        plt.scatter(x1[i], y1[i], color="red")

        # Perturbations
        # for j in range(1):
        #     plt.plot(x2s[j][:i + 1], y2s[j][:i + 1], color=colors[j])
        #     plt.scatter(x2s[j][i], y2s[j][i], color=colors[j])

        # Sun
        # plt.scatter(x2, y2, color="yellow")

        # H2 probes
        plt.plot(h2v1x[:i + 1], h2v1y[:i + 1], color="pink")
        plt.scatter(h2v1x[i], h2v1y[i], color="pink")

        plt.plot(h2v2x[:i + 1], h2v2y[:i + 1], color="orange")
        plt.scatter(h2v2x[i], h2v2y[i], color="orange")

        writer.grab_frame()
        plt.cla()
