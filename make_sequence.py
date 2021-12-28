
import csv
import colorsys
import time
import math
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from scipy.spatial.transform import Rotation as R
import numpy as np
import random

def hsl_to_rgb(h, s, l):
    h = h/360
    rgb = colorsys.hls_to_rgb(h, l, s)
    return [int(x*255) for x in rgb]

def getScale(coords):
    xs, ys, zs = zip(*coords)
    return min(max(xs)-min(xs), max(ys)-min(ys))

def centrifyCoords(coords):
    xs, ys, zs = zip(*coords)
    average_x = (min(xs)+max(xs))/2
    average_y = (min(ys)+max(ys))/2
    average_z = (min(zs)+max(zs))/2
    return [(x-average_x, y-average_y, z-average_z) for x, y, z in coords]

def distSqr(a, b):
    return math.pow(a[0]-b[0], 2) + math.pow(a[1]-b[1], 2) + math.pow(a[2]-b[2], 2)

def writeSequence(name, folder, sequence):
    num_lights = len(sequence[0])//3
    filename = "{}/night_{}.csv".format(folder, name)
    with open(filename, "w") as f:
        headings = ["FRAME_ID"]
        for i in range(num_lights):
            headings += ["{}_{}".format(char, i) for char in "RGB"]
        f.write(",".join(headings)+"\n")
        for i, frame in enumerate(sequence):
            row = [i] + frame
            f.write(",".join([str(x) for x in row])+"\n")

# coordinates and sequence folder for matt parker's xmas tree
coords_mattparker = []
with open("coords_2021.csv", encoding="utf-8-sig") as f:
    for row in csv.reader(f):
        row = [float(x) for x in row]
        coords_mattparker.append(row)
scale_mattparker = getScale(coords_mattparker)
folder_mattparker = "sequences"

# coordinates and sequence folder for harvard's xmas tree
coords_harvard = []
with open("coords_adjusted.txt") as f:
    for row in f.readlines():
        row = row.strip().strip("[]").split(", ")
        row = [float(x) for x in row]
        coords_harvard.append(row)
scale_harvard = getScale(coords_harvard)
folder_harvard = "sequences_harvard"

def cubeGridRotate(args):
    n, coords, scale, num_frames = args
    rainbow_size = scale*10.0
    colour_speed = 360/num_frames*3
    rotate_speed = 360/num_frames*2

    # rotation
    r = R.from_euler("xyz",
        [rotate_speed*n,
         rotate_speed*n,
         rotate_speed*n,
        ],
        degrees=True)
    new_coords = r.apply(coords)

    num_colours = 3
    grid_size = scale * (0.2 + 0.3*abs(math.sin(math.radians(rotate_speed*n))))
    assert 0.55 > grid_size/scale > 0.15

    frame = []
    for x, y, z in new_coords:
        colour = int(((x // grid_size) + (y // grid_size) + (z // grid_size)) % num_colours)
        if colour == 0:
            fade = (x%rainbow_size)/rainbow_size*360
        elif colour == 1:
            fade = (y%rainbow_size)/rainbow_size*360
        elif colour == 2:
            fade = (z%rainbow_size)/rainbow_size*360
        h = (fade + colour_speed*n + 60*colour) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def orbit(args):
    n, coords, scale, num_frames = args
    rainbow_size = scale*10.0
    colour_speed = 360/num_frames*3
    rotate_speed = 360/num_frames*2

    # rotation
    r = R.from_euler("xyz",
        [rotate_speed*n,
         rotate_speed*n,
         rotate_speed*n,
        ],
        degrees=True)
    new_coords = r.apply(coords)

    planets = [
        # position,       radius squared
        ((scale/2, 0, 0), math.pow(scale/5,2)),
        ((-scale/2, 0, 0), math.pow(scale/3,2)),
    ]

    #threshold = math.pow(scale/2, 2)

    num_colours = 3
    #factors = primeFactors(num_frames, num_colours)
    #colours = [(colour_speed*n*f + c*360/num_colours)%360 for c, f in enumerate(factors)]
    #grid_size = scale * (0.2 + 0.3*abs(math.sin(math.radians(rotate_speed*n))))
    #assert 0.55 > grid_size/scale > 0.15

    frame = []
    for x, y, z in new_coords:
        #colour = int(((x // grid_size) + (y // grid_size) + (z // grid_size)) % num_colours)

        for i, planet in enumerate(planets):
            planet_coords, planet_size = planet
            dist = distSqr((x, y, z), planet_coords)
            if dist <= planet_size:
                colour = i
                break
        else:
            frame.extend([0, 0, 0])
            continue

        #colour = random.randint(0, 2)

        if colour == 0:
            fade = (x%rainbow_size)/rainbow_size*360
        elif colour == 1:
            fade = (y%rainbow_size)/rainbow_size*360
        elif colour == 2:
            fade = (z%rainbow_size)/rainbow_size*360
        h = (fade + colour_speed*n + 60*colour) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def makeSequence(frameFunc, coords, scale, num_frames):
    coords = centrifyCoords(coords)
    start = timer()
    print("starting computations on {} cores".format(cpu_count()))
    values = [(n, coords, scale, num_frames) for n in range(num_frames)]

    with Pool() as pool:
        sequence = pool.map(frameFunc, values)
        print(len(sequence), type(sequence))
    end = timer()
    print("elapsed time: {}".format(end-start))

    """sequence = []
    for n in range(num_frames):
        frame = makeRainbowFrame((n, coords, scale, num_frames))
        sequence.append(frame)"""
    return sequence

num_frames = 120*6
frame_functions = [cubeGridRotate, orbit]
for frameFunc in frame_functions:
    seq_mattparker = makeSequence(frameFunc, coords_mattparker, scale_mattparker, num_frames)
    writeSequence(frameFunc.__name__, folder_mattparker, seq_mattparker)

    seq_harvard    = makeSequence(frameFunc, coords_harvard,    scale_harvard,    num_frames)
    writeSequence(frameFunc.__name__, folder_harvard, seq_harvard)

