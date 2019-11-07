
import json
import numpy as np
import matplotlib.pyplot as plt




def main():
    jsonFile = "../DataSet/PreprocessedData/intel_corrected_log"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        map = input['map']

    numSamplesPerRev = len(map[list(map)[0]]['range'])  # Get how many points per revolution
    angularStep = np.pi * 2 / numSamplesPerRev
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    xx = []
    yy = []
    for key in sorted(map.keys()):
        count += 1
        x, y, theta, range = map[key]['x'], map[key]['y'], map[key]['theta'], map[key]['range']
        rads = np.linspace(theta, theta + np.pi , num=numSamplesPerRev)
        range = np.asarray(range)
        range_idx = range < 50
        range = range[range_idx]
        rads = rads[range_idx]
        px = x + np.cos(rads) * range
        py = y + np.sin(rads) * range
        xx.append(x)
        yy.append(y)
        if count % 1 == 0:
            #plt.scatter(px, py, c= 'k', s = 1)

            plt.scatter(x, y, c='r', s=35)
        if count % 3 == 0:
            plt.plot(xx, yy)
            plt.show()

    plt.show()


if __name__ == '__main__':
    main()