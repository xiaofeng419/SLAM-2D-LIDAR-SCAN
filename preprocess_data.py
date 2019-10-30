import json

inputLogFile = "DataSet/intel.clf"
inputRelationFile = "DataSet/intel.relations"

map = {}
with open(inputLogFile, "r") as log_file:
    for line in log_file:
        if line.startswith('FLASER'):
            lineTokens = line.split()
            numPoints = int(lineTokens[1])
            range = lineTokens[2: numPoints + 2]
            range = [float(r) for r in range]
            x, y, theta, timeStamp = float(lineTokens[numPoints + 2]), float(lineTokens[numPoints + 3]), float(lineTokens[numPoints + 4]), float(lineTokens[numPoints + 8])
            map[timeStamp] = {'x': x, 'y': y, 'theta': theta, 'range': range}

relation_timeStamp1 = {}
relation_timeStamp2 = {}
with open(inputRelationFile, "r") as log_file:
    for line in log_file:
        lineTokens = line.split()
        x, y, theta, timeStamp1, timeStamp2 = float(lineTokens[2]), float(lineTokens[3]), float(lineTokens[7]), float(lineTokens[0]), float(lineTokens[1])
        relation_timeStamp1[timeStamp1] = {'x': x, 'y': y, 'theta': theta , 'timeStamp2': timeStamp2}
        relation_timeStamp2[timeStamp2] = {'x': x, 'y': y, 'theta': theta , 'timeStamp1': timeStamp1}

outputFile = inputLogFile[: inputLogFile.find('.')] + '_clean.json'
json_data = {
    'map': map,
    'relation_timeStamp1': relation_timeStamp1,
    'relation_timeStamp2': relation_timeStamp2,
}

with open(outputFile, 'w') as fp:
    json.dump(json_data, fp, sort_keys=True, indent=4)