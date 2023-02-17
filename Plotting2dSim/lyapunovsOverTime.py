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
f = "D:\\2dSunEarthChaos\\singleProbeData\\leOverTime.csv"

df = pd.read_csv(f)

h1y = df.loc[:, 'h1']
h1x = 0.001
h2y = df.loc[:, 'h2']
h2x = -0.001

colors = ['g', 'r', 'c', 'm', 'y', 'k', 'black', 'aqua', 'lime', 'maroon', 'navy', 'olive', 'purple', 'grey', 'silver',
          'teal', 'white', 'orange', 'hotpink', 'aquamarine']

with writer.saving(fig, "mp4s\\h2lyapunovTimeEvolution.mp4", 200):
    for i in range(len(h1y)):
        if i % 100 == 0:
            print(i)

        # h1
        # plt.scatter(h1x, h1y[i], color="blue")
        plt.scatter(h2x, h2y[i], color="red")

        writer.grab_frame()
        plt.cla()
