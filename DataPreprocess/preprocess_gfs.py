import json

inputLogFile = "../DataSet/RawData/csail.corrected.gfs"
outputFile = "../DataSet/PreprocessedData/csail_gfs"

map = {}
with open(inputLogFile, "r") as log_file:
    for line in log_file:
        if line.startswith('LASER_READING'):
            lineTokens = line.split()
            numPoints = int(lineTokens[1])
            range = lineTokens[2: numPoints + 2]
            range = [float(r) for r in range]
            x, y, theta, timeStamp = float(lineTokens[numPoints + 2]), float(lineTokens[numPoints + 3]), float(lineTokens[numPoints + 4]), float(lineTokens[numPoints + 5])
            map[timeStamp] = {'x': x, 'y': y, 'theta': theta, 'range': range}

json_data = {
    'map': map,
}

with open(outputFile, 'w') as fp:
    json.dump(json_data, fp, sort_keys=True, indent=4)