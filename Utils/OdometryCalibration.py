import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Utils.OccupancyGrid import OccupancyGrid
from scipy.ndimage import gaussian_filter
import math

def plotMove():
    a = 1
def main():
    jsonFile = "../DataSet/PreprocessedData/intel_gfs"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        sensorData = input['map']

    jsonFile = "../DataSet/PreprocessedData/intel_corrected_log"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        gtData = input['map']
    count = 0

    EstMoveList, GtMoveList, EstGtMoveList = [], [], []
    for key in sorted(sensorData.keys()):
        count += 1
        if count == 1:
            prevReading = sensorData[key]
            prevGtReading = gtData[key]
            continue
        reading = sensorData[key]
        gtReading = gtData[key]
        prevEstimatedX, prevEstimatedY, prevEstimatedTheta, prevRMeasure = prevReading['x'], prevReading['y'], prevReading['theta'], prevReading[
            'range']
        prevGtX, prevGtY, prevGtTheta, prevGtMeasure = prevGtReading['x'], prevGtReading['y'], prevGtReading['theta'], prevGtReading['range']
        estimatedX, estimatedY, estimatedTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        gtX, gtY, gtTheta, gtMeasure = gtReading['x'], gtReading['y'], gtReading['theta'], gtReading['range']
        dEstimatedMove = math.sqrt((estimatedX - prevEstimatedX) ** 2 + (estimatedY - prevEstimatedY) ** 2)
        dGtMove = math.sqrt((gtX - prevGtX) ** 2 + (gtY - prevGtY) ** 2)
        EstMoveList.append(dEstimatedMove)
        GtMoveList.append(dGtMove)
        EstGtMoveList.append(dEstimatedMove - dGtMove)
        prevGtReading = gtReading
        prevReading = reading
        print(count)

    EstGtMoveArray = np.asarray(EstGtMoveList)
    EstMoveSortedIdx = np.argsort(EstMoveList)
    EstMoveSortedArray = np.sort(EstMoveList)
    EstGtMoveSortedArray = EstGtMoveArray[EstMoveSortedIdx]
    plt.scatter(EstMoveSortedArray, EstGtMoveSortedArray)
    plt.show()
    GtMoveArray = np.asarray(GtMoveList)
    GtMoveSortedArray = GtMoveArray[EstMoveSortedIdx]
    plt.plot(EstMoveSortedArray)
    plt.plot(GtMoveSortedArray)
    plt.show()
if __name__ == '__main__':
    main()