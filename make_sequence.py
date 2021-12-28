
import csv
import colorsys
import time
import math
from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from scipy.spatial.transform import Rotation as R
import numpy as np

primes = [2, 3, 5, 7, 11]

mem = {}
def primeFactors(n, count):
    if (n, count) in mem:
        return mem[(n, count)]
    factors = []
    for prime in primes:
        if n % prime == 0:
            factors.append(prime)
        if len(factors) >= count:
            break
    if len(factors) < count:
        raise Exception("this isn't a very good frame count")
    mem[(n, count)] = factors
    return factors

def hsl_to_rgb(h, s, l):
    h = h/360
    rgb = colorsys.hls_to_rgb(h, l, s)
    return [int(x*255) for x in rgb]

def getScale(coords):
    xs, ys, zs = zip(*coords)
    return min(max(xs)-min(xs), max(ys)-min(ys))

def centrifyCoords(coords):
    xs, ys, zs = zip(*coords)
    average_x = (min(xs)+max(xs))/2#np.median(xs)
    average_y = (min(ys)+max(ys))/2#np.median(ys)
    average_z = (min(zs)+max(zs))/2#np.median(zs)
    return [(x-average_x, y-average_y, z-average_z) for x, y, z in coords]

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

def makeRainbowFrame(args):#n, coords, scale, num_frames):
    n, coords, scale, num_frames = args
    rainbow_size = scale*0.5
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

    #threshold = math.pow(scale/2, 2)

    num_colours = 2
    #factors = primeFactors(num_frames, num_colours)
    #colours = [(colour_speed*n*f + c*360/num_colours)%360 for c, f in enumerate(factors)]
    grid_size = scale * (0.3 + 0.2*math.sin(math.radians(rotate_speed*n)))
    assert 0.55 > grid_size/scale > 0.05

    frame = []
    for x, y, z in new_coords:
        #dist = math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2)
        #if dist > threshold:
        #    frame.extend([0, 0, 0])
        #    continue
        colour = int(((x // grid_size) + (y // grid_size) + (z // grid_size)) % num_colours)
        #h = colours[colour]
        #h = (90*colour) % 360
        h = (colour_speed*n + 60*colour) % 360
        #h = ((z%rainbow_size)/rainbow_size*360 + colour_speed*n + 180*extra) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def makeRainbow(coords, scale, num_frames):
    coords = centrifyCoords(coords)
    start = timer()
    print("starting computations on {} cores".format(cpu_count()))
    values = [(n, coords, scale, num_frames) for n in range(num_frames)]

    with Pool() as pool:
        sequence = pool.map(makeRainbowFrame, values)
        print(len(sequence), type(sequence))
    end = timer()
    print("elapsed time: {}".format(end-start))

    """sequence = []
    for n in range(num_frames):
        frame = makeRainbowFrame((n, coords, scale, num_frames))
        sequence.append(frame)"""
    return sequence

num_frames = 120*6
rainbow_mattparker = makeRainbow(coords_mattparker, scale_mattparker, num_frames)
rainbow_harvard    = makeRainbow(coords_harvard, scale_harvard, num_frames)

writeSequence("rainbow", folder_mattparker, rainbow_mattparker)
writeSequence("rainbow", folder_harvard, rainbow_harvard)

