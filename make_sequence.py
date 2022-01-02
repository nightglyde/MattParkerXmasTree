
import csv
import colorsys
import time
import math
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from scipy.spatial.transform import Rotation as R
import numpy as np
import random

# set max brightness because tree power supply struggles
MAX_BRIGHTNESS = 50

def hsl_to_rgb(h, s, l):
    h = h/360
    rgb = colorsys.hls_to_rgb(h, l, s)
    return [int(x*MAX_BRIGHTNESS) for x in rgb]

def combine(a, b):
    return [max(pair) for pair in zip(a, b)]

def combineFrames(a, b):
    num_lights = len(a)//3
    frame = []
    for i in range(num_lights):
        ra, ga, ba = a[i*3], a[i*3+1], a[i*3+2]
        rb, gb, bb = b[i*3], b[i*3+1], b[i*3+2]

        if ra or ga or ba:
            frame.extend([ra, ga, ba])
        else:
            frame.extend([rb, gb, bb])
    return frame

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
    colour_speed = 360/num_frames*2
    rotate_speed = 360/num_frames

    # rotation
    r = R.from_euler("xyz",
        [rotate_speed*n*2,
         rotate_speed*n*3,
         rotate_speed*n*5,
        ],
        degrees=True)
    new_coords = r.apply(coords)

    num_colours = 3
    grid_size = scale * 0.5#(0.3 + 0.2*abs(math.sin(math.radians(rotate_speed*n))))
    assert 0.55 > grid_size/scale > 0.25

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

def cubeGridRotatePlanes(args):
    n, coords, scale, num_frames = args
    rainbow_size = scale*5.0
    colour_speed = 360/num_frames*2
    rotate_speed = 360/num_frames

    # rotation
    r = R.from_euler("xyz",
        [rotate_speed*n*2,
         rotate_speed*n*3,
         rotate_speed*n*5,
        ],
        degrees=True)
    new_coords = r.apply(coords)

    grid_size = scale * 0.5#(0.3 + 0.2*abs(math.sin(math.radians(rotate_speed*n))))
    assert 0.55 > grid_size/scale > 0.25

    threshold = scale*0.05

    frame = []
    for x, y, z in new_coords:
        x_plane = x%grid_size < threshold or (grid_size-threshold) < x%grid_size
        y_plane = y%grid_size < threshold or (grid_size-threshold) < y%grid_size
        z_plane = z%grid_size < threshold or (grid_size-threshold) < z%grid_size

        if not (x_plane or y_plane or z_plane):
            frame.extend([0, 0, 0])
            continue

        x_h = ((x%rainbow_size)/rainbow_size*360 + colour_speed*n +   0) % 360
        y_h = ((y%rainbow_size)/rainbow_size*360 + colour_speed*n +  60) % 360
        z_h = ((z%rainbow_size)/rainbow_size*360 + colour_speed*n + 120) % 360
        s = 1.0
        l = 0.5
        x_colour = hsl_to_rgb(x_h, s, l)
        y_colour = hsl_to_rgb(y_h, s, l)
        z_colour = hsl_to_rgb(z_h, s, l)

        rgb = [0, 0, 0]
        if x_plane:
            rgb = combine(rgb, x_colour)
        if y_plane:
            rgb = combine(rgb, y_colour)
        if z_plane:
            rgb = combine(rgb, z_colour)
        frame.extend(rgb)
    return frame

def binaryStars(args):
    n, coords, scale, num_frames = args
    rainbow_size = scale*10.0
    colour_speed = 360/num_frames*2
    rotate_speed = 360/num_frames

    # rotation
    r = R.from_euler("xyz",
        [0,#rotate_speed*n*2,
         rotate_speed*n*3,
         rotate_speed*n*4,
        ],
        degrees=True)
    new_coords = r.apply(coords)

    spheres = [
        # position,       radius squared
        ((scale*0.5, 0, 0), math.pow(scale*0.55,2)),
        ((scale*-0.6655, 0, 0), math.pow(scale*0.5,2)),
    ]

    frame = []
    for x, y, z in new_coords:
        for i, sphere in enumerate(spheres):
            sphere_coords, sphere_size = sphere
            dist = distSqr((x, y, z), sphere_coords)
            if dist <= sphere_size:
                colour = i
                break
        else:
            frame.extend([0, 0, 0])
            continue

        fade = (x%rainbow_size)/rainbow_size*360
        h = (fade + colour_speed*n + 60*colour) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def translateZ(coords, offset):
    new_coords = []
    for x, y, z in coords:
        new_coords.append((x, y, z+offset))
    return new_coords

def spherical(args):
    n, coords, scale, num_frames = args
    colour_speed = 360/num_frames*10
    rotate_speed = 360/num_frames

    # rotation
    r = R.from_euler("xyz",
        [0,#rotate_speed*n,
         0,#rotate_speed*n,
         rotate_speed*n*3,
        ],
        degrees=True)

    offset = scale * 0.5 * math.sin(math.radians(rotate_speed*n))

    new_coords = translateZ(r.apply(coords), offset)

    spheres = [
        # position,       radius squared
        ((scale*0.2, 0, 0), math.pow(scale*0.2,2)),
    ]

    frame = []
    for x, y, z in new_coords:
        for i, sphere in enumerate(spheres):
            sphere_coords, sphere_size = sphere
            dist = distSqr((x, y, z), sphere_coords)
            if dist <= sphere_size:
                colour = i
                break
        else:
            frame.extend([0, 0, 0])
            continue

        #fade = (x%rainbow_size)/rainbow_size*360
        h = (colour_speed*n + 60*colour) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def makeSequence(frameFunc, coords, scale, num_frames):
    coords = centrifyCoords(coords)
    start = timer()
    print("generating sequence {} on {} cores".format(frameFunc.__name__, cpu_count()))
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

def snakeyFrame(args):
    n, coords, orig_sequence, snakey_frames = args

    num_frames = len(orig_sequence)

    frame = [0 for i in range(len(coords)*3)]
    for s in range(snakey_frames):
        f = (n+s)%num_frames
        frame = combineFrames(orig_sequence[f], frame)
    return frame

def snakey(frameFunc, coords, scale, num_frames, snakey_frames):
    coords = centrifyCoords(coords)
    start = timer()
    print("generating snakey sequence {} on {} cores".format(frameFunc.__name__, cpu_count()))
    values = [(n, coords, scale, num_frames) for n in range(num_frames)]

    with Pool() as pool:
        orig_sequence = pool.map(frameFunc, values)
        print(len(orig_sequence), type(orig_sequence))
    end = timer()
    print("elapsed time: {}".format(end-start))

    start = timer()
    print("snakifying sequence on {} cores".format(cpu_count()))
    values2 = [(n, coords, orig_sequence, snakey_frames) for n in range(num_frames)]

    with Pool() as pool:
        sequence = pool.map(snakeyFrame, values2)
        print(len(orig_sequence), type(orig_sequence))
    end = timer()
    print("elapsed time: {}".format(end-start))

    return sequence

num_frames = 120*6
frame_functions = [
    cubeGridRotate,
    cubeGridRotatePlanes,
    binaryStars,
]
for frameFunc in frame_functions:
    seq_mattparker = makeSequence(frameFunc, coords_mattparker, scale_mattparker, num_frames)
    writeSequence(frameFunc.__name__, folder_mattparker, seq_mattparker)

    seq_harvard    = makeSequence(frameFunc, coords_harvard,    scale_harvard,    num_frames)
    writeSequence(frameFunc.__name__, folder_harvard, seq_harvard)

snakey_frames = num_frames//6
snakey_frame_functions = [
    spherical,
]
for frameFunc in snakey_frame_functions:
    seq_mattparker = snakey(frameFunc, coords_mattparker, scale_mattparker, num_frames, snakey_frames)
    writeSequence("snakey"+frameFunc.__name__, folder_mattparker, seq_mattparker)

    seq_harvard    = snakey(frameFunc, coords_harvard, scale_harvard, num_frames, snakey_frames)
    writeSequence("snakey"+frameFunc.__name__, folder_harvard, seq_harvard)

