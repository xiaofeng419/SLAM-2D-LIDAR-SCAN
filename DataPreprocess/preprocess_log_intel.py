import json

refTimeJsonFile = "../DataSet/PreprocessedData/intel_raw_refTime"
with open(refTimeJsonFile, 'r') as f:
    input = json.load(f)
    mapRef = input['map']

mapRefrange = {}
for key in mapRef:
    #mapRefxytheta[(mapRef[key]['x'], mapRef[key]['y'], mapRef[key]['theta'])] = key
    mapRefrange[tuple(mapRef[key]['range'])] = key

inputGfsFile = "../DataSet/RawData/intel.gfs"
inputCorrectedLogFile = "../DataSet/RawData/intel.gfs.log"
outputGfsFile = "../DataSet/PreprocessedData/intel_gfs"
outputCorrectedLogFile = "../DataSet/PreprocessedData/intel_corrected_log"

mapGfs = {}  # non-corrected
mapLog = {}  # corrected
gfsTime2RefTimeMap = {}

with open(inputGfsFile, "r") as gfs_file:
    for line in gfs_file:
        if line.startswith('LASER_READING'):
            lineTokens = line.split()
            numPoints = int(lineTokens[1])
            range = lineTokens[2: numPoints + 2]
            range = [float(r) for r in range]
            x, y, theta, gfsTimeStamp = float(lineTokens[numPoints + 2]), float(lineTokens[numPoints + 3]), float(lineTokens[numPoints + 4]), float(lineTokens[numPoints + 5])
            refTimeStamp = mapRefrange[tuple(range)]
            gfsTime2RefTimeMap[gfsTimeStamp] = refTimeStamp
            mapGfs[refTimeStamp] = {'x': x, 'y': y, 'theta': theta, 'range': range}
json_gfs_data = {
    'map': mapGfs,
}
with open(outputGfsFile, 'w') as fp:
    json.dump(json_gfs_data, fp, sort_keys=True, indent=4)

with open(inputCorrectedLogFile, "r") as log_file:
    for line in log_file:
        if line.startswith('FLASER'):
            lineTokens = line.split()
            numPoints = int(lineTokens[1])
            range = lineTokens[2: numPoints + 2]
            range = [float(r) for r in range]
            x, y, theta, gfsTimeStamp = float(lineTokens[numPoints + 2]), float(lineTokens[numPoints + 3]), float(lineTokens[numPoints + 4]), float(lineTokens[numPoints + 8])
            refTimeStamp = gfsTime2RefTimeMap[gfsTimeStamp]
            mapLog[refTimeStamp] = {'x': x, 'y': y, 'theta': theta, 'range': range}

json_log_data = {
    'map': mapLog,
}

with open(outputCorrectedLogFile, 'w') as fp:
    json.dump(json_log_data, fp, sort_keys=True, indent=4)