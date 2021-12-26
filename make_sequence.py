
import csv
import colorsys

def hsl_to_rgb(h, s, l):
    h = h/360
    rgb = colorsys.hls_to_rgb(h, l, s)
    return [int(x*255) for x in rgb]

def getScale(coords):
    x, y, z = coords[0]
    xmin = xmax = x
    ymin = ymax = y
    for x, y, z in coords:
        xmin = min(x, xmin)
        xmax = max(x, xmax)
        ymin = min(y, ymin)
        ymax = max(y, ymax)
    rangex = xmax - xmin
    rangey = ymax - ymin
    return min(rangex, rangey)

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

def makeRainbowFrame(n, coords, scale, num_frames):
    rainbow_size = scale*2.0
    colour_change_rate = 360/num_frames
    frame = []
    for x, y, z in coords:
        h = ((z%rainbow_size)/rainbow_size*360 + colour_change_rate*n) % 360
        s = 1.0
        l = 0.5
        rgb = hsl_to_rgb(h, s, l)
        frame.extend(rgb)
    return frame

def makeRainbow(coords, scale, num_frames):
    sequence = []
    for n in range(num_frames):
        frame = makeRainbowFrame(n, coords, scale, num_frames)
        sequence.append(frame)
    return sequence

num_frames = 120
rainbow_mattparker = makeRainbow(coords_mattparker, scale_mattparker, num_frames)
rainbow_harvard    = makeRainbow(coords_harvard, scale_harvard, num_frames)

writeSequence("rainbow", folder_mattparker, rainbow_mattparker)
writeSequence("rainbow", folder_harvard, rainbow_harvard)

