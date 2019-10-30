import json


inputRelationFile = "DataSet/intel.relations"
outputFile = "DataSet/intel_relation_processed"

relation_timeStamp1 = {}
relation_timeStamp2 = {}
with open(inputRelationFile, "r") as log_file:
    for line in log_file:
        lineTokens = line.split()
        x, y, theta, timeStamp1, timeStamp2 = float(lineTokens[2]), float(lineTokens[3]), float(lineTokens[7]), float \
            (lineTokens[0]), float(lineTokens[1])
        relation_timeStamp1[timeStamp1] = {'x': x, 'y': y, 'theta': theta , 'timeStamp2': timeStamp2}
        relation_timeStamp2[timeStamp2] = {'x': x, 'y': y, 'theta': theta , 'timeStamp1': timeStamp1}

json_data = {
    'relation_timeStamp1': relation_timeStamp1,
    'relation_timeStamp2': relation_timeStamp2,
}

with open(outputFile, 'w') as fp:
    json.dump(json_data, fp, sort_keys=True, indent=4)