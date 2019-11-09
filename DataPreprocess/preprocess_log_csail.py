import json

refTimeJsonFile = "../DataSet/PreprocessedData/csail_gfs"
with open(refTimeJsonFile, 'r') as f:
    input = json.load(f)
    mapRef = input['map']

refTimeList = sorted(mapRef.keys())

inputLogFile = "../DataSet/RawData/csail.corrected.log"
outputFile = "../DataSet/PreprocessedData/csail_corrected_log"

map = {}
count = 0
with open(inputLogFile, "r") as log_file:
    for line in log_file:
        if line.startswith('FLASER'):
            lineTokens = line.split()
            numPoints = int(lineTokens[1])
            range = lineTokens[2: numPoints + 2]
            range = [float(r) for r in range]
            x, y, theta, timeStamp = float(lineTokens[numPoints + 2]), float(lineTokens[numPoints + 3]), float(lineTokens[numPoints + 4]), float(lineTokens[numPoints + 8])
            timeStamp = refTimeList[count]
            count += 1
            map[timeStamp] = {'x': x, 'y': y, 'theta': theta, 'range': range}

json_data = {
    'map': map,
}

with open(outputFile, 'w') as fp:
    json.dump(json_data, fp, sort_keys=True, indent=4)