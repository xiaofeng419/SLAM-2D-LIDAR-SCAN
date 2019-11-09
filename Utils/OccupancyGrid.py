import json
import numpy as np
import matplotlib.pyplot as plt

def spokesGrid(radius, unitGridSize):
    numGrids = 2 * radius



def main():
    jsonFile = "../DataSet/PreprocessedData/intel_corrected_log"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        map = input['map']

    numSamplesPerRev = len(map[list(map)[0]]['range'])  # Get how many points per revolution
    angularStep = np.pi / numSamplesPerRev
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    xx = []
    yy = []
    for key in sorted(map.keys()):
        count += 1
        x, y, theta, range = map[key]['x'], map[key]['y'], map[key]['theta'], map[key]['range']
        rads = np.linspace(theta - np.pi / 2, theta + np.pi / 2 , num=numSamplesPerRev)
        range = np.asarray(range)
        range_idx = range < 50
        range = range[range_idx]
        rads = rads[range_idx]
        px = x + np.cos(rads) * range
        py = y + np.sin(rads) * range
        xx.append(x)
        yy.append(y)
        if count % 1 == 0:
            if count == 1:
                startx, starty = x, y
            plt.scatter(x, y, c='r', s=35)
            plt.scatter(px, py, c='k', s=1)
    plt.scatter(startx, starty, c='g', s=500)
    plt.scatter(x, y, c='b', s=500)
    plt.plot(xx, yy)
    plt.show()


if __name__ == '__main__':
    main()