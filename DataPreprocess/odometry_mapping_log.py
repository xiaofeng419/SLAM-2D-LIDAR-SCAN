import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def updateTrajectoryPlot(matchedReading, xTrajectory, yTrajectory, colors, count):
    x, y, theta, range = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
    xTrajectory.append(x)
    yTrajectory.append(y)
    if count % 1 == 0:
        plt.scatter(x, y, color=next(colors), s=35)

def main():
    jsonFile = "../DataSet/PreprocessedData/intel_gfs"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        map = input['map']

    numSamplesPerRev = len(map[list(map)[0]]['range'])  # Get how many points per revolution
    angularStep = np.pi / numSamplesPerRev
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    xx = []
    yy = []
    colors = iter(cm.rainbow(np.linspace(1, 0, len(map) + 1)))
    for key in sorted(map.keys()):
        count += 1
        x, y, theta, range = map[key]['x'], map[key]['y'], map[key]['theta'], map[key]['range']
        rads = np.linspace(theta - np.pi / 2, theta + np.pi / 2, num=numSamplesPerRev)
        range = np.asarray(range)
        range_idx = range < 10
        range = range[range_idx]
        rads = rads[range_idx]
        px = x + np.cos(rads) * range
        py = y + np.sin(rads) * range
        xx.append(x)
        yy.append(y)
        plt.scatter(px, py, color='black', s=1)
        if count == 1:
            startx, starty = x, y
        updateTrajectoryPlot(map[key], xx, yy, colors, count)

    plt.scatter(startx, starty, c='r', s=500)
    plt.scatter(x, y, c='b', s=500)
    plt.plot(xx, yy)
    plt.show()


if __name__ == '__main__':
    main()