import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Utils.OccupancyGrid import OccupancyGrid
from scipy.ndimage import gaussian_filter
import math

def plotMove(rawMinusGtMoveList, rawMoveList, gtMoveList):
    EstGtMoveArray = np.asarray(rawMinusGtMoveList)
    EstMoveSortedIdx = np.argsort(rawMoveList)
    EstMoveSortedArray = np.sort(rawMoveList)
    EstGtMoveSortedArray = EstGtMoveArray[EstMoveSortedIdx]
    plt.scatter(EstMoveSortedArray, EstGtMoveSortedArray)
    plt.xlabel("m")
    plt.ylabel("m")
    plt.show()
    GtMoveArray = np.asarray(gtMoveList)
    GtMoveSortedArray = GtMoveArray[EstMoveSortedIdx]
    plt.plot(EstMoveSortedArray)
    plt.plot(GtMoveSortedArray)
    plt.xlabel("m")
    plt.ylabel("m")
    plt.show()

def plotTheta(errorTheta1List, rawMoveList):
    rawMoveSortedIdx = np.argsort(rawMoveList)
    rawMoveSortedArray = np.sort(rawMoveList)
    errorTheta1Array = np.asarray(errorTheta1List)
    errorSortedTheta1Array = errorTheta1Array[rawMoveSortedIdx]
    plt.scatter(rawMoveSortedArray, errorSortedTheta1Array, s = 1)
    plt.xlabel("m")
    plt.ylabel("deg")
    plt.show()



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

    rawMoveList, gtMoveList, rawMinusGtMoveList, turningAngleList = [], [], [], []
    errorTheta1List, errorTheta2List, errorTheta2AsDirList = [], [], []
    for key in sorted(sensorData.keys()):
        count += 1
        if count == 1:
            prevReading = sensorData[key]
            prevGtReading = gtData[key]
            prevRawThetaM = None
            continue
        reading = sensorData[key]
        gtReading = gtData[key]
        prevRawX, prevRawY, prevRawTheta, prevRawRMeasure = prevReading['x'], prevReading['y'], prevReading['theta'], prevReading[
            'range']
        prevGtX, prevGtY, prevGtTheta, prevGtMeasure = prevGtReading['x'], prevGtReading['y'], prevGtReading['theta'], prevGtReading['range']
        rawX, rawY, rawTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        gtX, gtY, gtTheta, gtMeasure = gtReading['x'], gtReading['y'], gtReading['theta'], gtReading['range']

        # get move
        rawXMove, rawYMove, gtXMove, gtYMove = rawX - prevRawX, rawY - prevRawY, gtX - prevGtX, gtY - prevGtY
        rawMove = math.sqrt((rawX - prevRawX) ** 2 + (rawY - prevRawY) ** 2)
        gtMove = math.sqrt((gtX - prevGtX) ** 2 + (gtY - prevGtY) ** 2)

        if abs(rawMove) < 0.01:
            continue

        # get theta 1
        ## raw thetaM
        if rawMove > 0.3:
            if prevRawThetaM != None:

                if rawYMove > 0:
                    rawThetaM = math.acos(rawXMove / rawMove)
                else:
                    rawThetaM = -math.acos(rawXMove / rawMove)
                rawTheta1 = rawThetaM - prevRawThetaM
                if rawTheta1 < 0:
                    rawTheta1 = 2 * np.pi + rawTheta1
                if gtYMove > 0:
                    gtThetaM = math.acos(gtXMove / gtMove)
                else:
                    gtThetaM = - math.acos(gtXMove / gtMove)
                gtTheta1 = gtThetaM - prevGtThetaM
                if gtTheta1 < 0:
                    gtTheta1 = 2 * np.pi + gtTheta1

                errorTheta1 = rawTheta1 - gtTheta1
                if errorTheta1 > np.pi:
                   errorTheta1 = errorTheta1 - np.pi * 2
                elif errorTheta1 < - np.pi:
                   errorTheta1 = errorTheta1 + np.pi * 2

                prevRawThetaM = rawThetaM
                prevGtThetaM = gtThetaM
            else:
                if rawYMove > 0:
                    rawThetaM = math.acos(rawXMove / rawMove)
                else:
                    rawThetaM = -math.acos(rawXMove / rawMove)
                if gtYMove > 0:
                    gtThetaM = math.acos(gtXMove / gtMove)
                else:
                    gtThetaM = - math.acos(gtXMove / gtMove)
                prevRawThetaM = rawThetaM
                prevGtThetaM = gtThetaM
                errorTheta1 = None
                rawTheta1 = None
        else:
            prevRawThetaM = None
            prevGtThetaM = None
            errorTheta1 = None
            rawTheta1 = None


        # theta 2
        rawTheta2 = rawTheta - prevRawTheta
        gtTheta2 = gtTheta - prevGtTheta
        errorTheta2 = rawTheta2 - gtTheta2
        if errorTheta2 > np.pi:
            errorTheta2 = errorTheta2 - np.pi * 2
        elif errorTheta2 < - np.pi:
            errorTheta2 = errorTheta2 + np.pi * 2

        # use theta2 as moving direction
        #errorTheta2AsDir = rawTheta - gtThetaM



        # List appending
        rawMoveList.append(rawMove)
        gtMoveList.append(gtMove)
        rawMinusGtMoveList.append(rawMove - gtMove)
        errorTheta1List.append(errorTheta1)
        errorTheta2List.append(errorTheta2)

        turningAngleList.append(rawTheta1)
        #errorTheta2AsDirList.append(errorTheta2AsDir)

        print(count)
        prevGtReading = gtReading
        prevReading = reading

    rawMoveArray = np.asarray(rawMoveList)
    turningAngleArray = np.asarray(turningAngleList)
    turningAngleArray = turningAngleArray[rawMoveArray > 0.1]
    errorTheta1Array = np.asarray(errorTheta1List)
    errorTheta1Array = errorTheta1Array[rawMoveArray > 0.1]
    plt.scatter(turningAngleArray, errorTheta1Array)
    plt.show()

    plotMove(rawMinusGtMoveList, rawMoveList, gtMoveList)
    plotTheta(errorTheta1List, rawMoveList)
    plotTheta(errorTheta2List, rawMoveList)
    #plotTheta(errorTheta2AsDirList, rawMoveList)

if __name__ == '__main__':
    main()