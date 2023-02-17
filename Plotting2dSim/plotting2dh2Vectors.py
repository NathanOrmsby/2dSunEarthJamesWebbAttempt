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
writer = FFMpegWriter(fps=60, metadata=metadata)

fig = plt.figure()
l, = plt.plot([], [], 'k-')

# Read from csv files
# Desktop
f0 = "D:\\2dSunEarthChaos\\singleProbeData\\NormalizedGramSchmidt.csv"
f1 = "D:\\2dSunEarthChaos\\singleProbeData\\gramSchmidtVectors.csv"


df0 = pd.read_csv(f0)
df1 = pd.read_csv(f1)

# h2 Normalized:
u1x = df0.loc[:, "u1x"]
u1y = df0.loc[:, "u1y"]
u2x = df0.loc[:, "u2x"]
u2y = df0.loc[:, "u2y"]

# h2 Not normalized
v1x = df1.loc[:, "v1x"]
v1y = df1.loc[:, "v1y"]
v2x = df1.loc[:, "v2x"]
v2y = df1.loc[:, "v2y"]


# vector1 = np.array([[0, 0, u1x[i], u1y[i]] for i in range(len(u1x))])
# X, Y, U, V = zip(*vector1)
#
# vector2 = np.array([[0, 0, u2x[i], u2y[i]] for i in range(len(u2x))])
# x, y, u, v = zip(*vector2)

vector1 = np.array([[0, 0, v1x[i], v1y[i]] for i in range(len(v1x))])
X, Y, U, V = zip(*vector1)

vector2 = np.array([[0, 0, v2x[i], v2y[i]] for i in range(len(v2x))])
x, y, u, v = zip(*vector2)



colors = ['g', 'r', 'c', 'm', 'y', 'k', 'black', 'aqua', 'lime', 'maroon', 'navy', 'olive', 'purple', 'grey', 'silver', 'teal', 'white', 'orange', 'hotpink', 'aquamarine']

with writer.saving(fig, "mp4s\\gramSchmidtEvolution.mp4", 200):
    for i in range(len(u1x)):
        if i % 100 == 0:
            print(i)

        plt.quiver(X[i], Y[i], U[i], V[i], color='r')
        plt.quiver(x[i], y[i], u[i], v[i], color='g')

        writer.grab_frame()
        plt.cla()
