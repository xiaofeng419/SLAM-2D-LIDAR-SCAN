import json
import numpy as np
import matplotlib.pyplot as plt


jsonFile = "DataSet/intel_clean.json"

with open(jsonFile, 'r') as f:
    input = json.load(f)
    map = input['map']
    relation_timeStamp1 = input['relation_timeStamp1']
    relation_timeStamp2 = input['relation_timeStamp2']

numSamplesPerRev = len(map[list(map)[0]]['range'])  # Get how many points per revolution
angularStep = np.pi * 2 / numSamplesPerRev
count = 0
plt.figure(figsize=(19.20, 10.80))

for key in sorted(map.keys()):
    count += 1
    x, y, theta, range = map[key]['x'], map[key]['y'], map[key]['theta'], map[key]['range']
    rads = np.linspace(theta, theta + np.pi * 2, num=numSamplesPerRev)
    range = np.asarray(range)
    px =   np.cos(rads) * range
    py =   np.sin(rads) * range
    if count % 1 == 0 :
        #plt.scatter(px, py, c= 'k', s = 1)
        plt.scatter(x, y, c = 'r', s = 1)
plt.show()


