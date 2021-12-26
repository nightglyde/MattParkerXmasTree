
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv, os

print("""
Remember that the animation on the real tree might be significantly
faster or slower than the animation speed of this program.
""")

# use harvard's tree or mattparker's tree
# default is False
USE_HARVARD = False

# coordinates and sequence folder for matt parker's xmas tree
coords_mattparker = []
with open("coords_2021.csv", encoding="utf-8-sig") as f:
    for row in csv.reader(f):
        row = [float(x) for x in row]
        coords_mattparker.append(row)
folder_mattparker = "sequences"

# coordinates and sequence folder for harvard's xmas tree
coords_harvard = []
with open("coords_adjusted.txt") as f:
    for row in f.readlines():
        row = row.strip().strip("[]").split(", ")
        row = [float(x) for x in row]
        coords_harvard.append(row)
folder_harvard = "sequences_harvard"

# extract sequence from file
def get_sequence(filename):
    sequence_data = []
    with open(filename) as f:
        for row in csv.reader(f):
            sequence_data.append(row)

    frames = []
    for i, frame_data in enumerate(sequence_data[1:]): # skip heading row
        frame_num = int(frame_data[0])
        assert frame_num == i
        frame_data = [float(x)/255 for x in frame_data[1:]] # skip frame num
        frame = []
        for i in range(0, len(frame_data), 3):
            frame.append(frame_data[i:i+3])
        frames.append(frame)
    return frames

if USE_HARVARD:
    coords = coords_harvard
    folder = folder_harvard
else:
    coords  = coords_mattparker
    folder = folder_mattparker

sequences = []
for filename in os.listdir(folder):
    filename = "{}/{}".format(folder, filename)
    sequences.append((get_sequence(filename), filename))

# handle coords of lights on xmas tree
xs = []
ys = []
zs = []
max_coords = [0.0, 0.0, 0.0]
min_coords = [0.0, 0.0, 0.0]
for row in coords:
        x, y, z = row
        xs.append(x)
        ys.append(y)
        zs.append(z)
        for i, coord in enumerate(row):
            min_coords[i] = min(coord, min_coords[i])
            max_coords[i] = max(coord, max_coords[i])

ranges = [maxi - mini for mini, maxi in zip(min_coords, max_coords)]

# make the plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_box_aspect(ranges)

def frame_generator(seqs, repeats=1):
    while True:
        for seq, name in seqs:
            print(name)
            for repeat in range(repeats):
                for frame in seq:
                    yield frame

framegen = frame_generator(sequences, 3)

def update_colours(frame):
    colour_values = next(framegen)
    ax.clear()
    ax.scatter(xs, ys, zs, c=colour_values)

ani = animation.FuncAnimation(fig, update_colours, interval=17)

plt.show()

