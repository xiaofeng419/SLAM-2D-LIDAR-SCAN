import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Utils.OccupancyGrid import OccupancyGrid
from scipy.ndimage import gaussian_filter
import math
class ScanMatcher:
    def __init__(self, og, searchRadius, searchHalfRad, scanSigmaInNumGrid, moveRSigma, maxMoveDeviation, turnSigma, missMatchProbAtCoarse, coarseFactor):
        self.searchRadius = searchRadius
        self.searchHalfRad = searchHalfRad
        self.og = og
        self.scanSigmaInNumGrid = scanSigmaInNumGrid
        self.coarseFactor = coarseFactor
        self.moveRSigma = moveRSigma
        self.turnSigma = turnSigma
        self.missMatchProbAtCoarse = missMatchProbAtCoarse
        self.maxMoveDeviation = maxMoveDeviation

    def frameSearchSpace(self, estimatedX, estimatedY, unitLength, sigma, missMatchProbAtCoarse):
        maxScanRadius = 1.1 * self.og.lidarMaxRange + self.searchRadius
        xRangeList = [estimatedX - maxScanRadius, estimatedX + maxScanRadius]
        yRangeList = [estimatedY - maxScanRadius, estimatedY + maxScanRadius]
        idxEndX, idxEndY = int((xRangeList[1] - xRangeList[0]) / unitLength),  int((yRangeList[1] - yRangeList[0]) / unitLength)
        searchSpace = math.log(missMatchProbAtCoarse) * np.ones((idxEndY + 1, idxEndX + 1))

        self.og.checkAndExapndOG(xRangeList, yRangeList)
        xRangeListIdx, yRangeListIdx = self.og.convertRealXYToMapIdx(xRangeList, yRangeList)
        ogMap = self.og.occupancyGridVisited[yRangeListIdx[0]: yRangeListIdx[1], xRangeListIdx[0]: xRangeListIdx[1]] /\
                      self.og.occupancyGridTotal[yRangeListIdx[0]: yRangeListIdx[1], xRangeListIdx[0]: xRangeListIdx[1]]
        ogMap = ogMap > 0.5
        ogX = self.og.OccupancyGridX[yRangeListIdx[0]: yRangeListIdx[1], xRangeListIdx[0]: xRangeListIdx[1]]
        ogY = self.og.OccupancyGridY[yRangeListIdx[0]: yRangeListIdx[1], xRangeListIdx[0]: xRangeListIdx[1]]

        ogX, ogY = ogX[ogMap], ogY[ogMap]
        ogIdx = self.convertXYToSearchSpaceIdx(ogX, ogY, xRangeList[0], yRangeList[0], unitLength)
        searchSpace[ogIdx[1], ogIdx[0]] = 0
        probSP = self.generateProbSearchSpace(searchSpace, sigma)
        return xRangeList, yRangeList, probSP

    def generateProbSearchSpace(self, searchSpace, sigma):
        probSP = gaussian_filter(searchSpace, sigma=sigma)
        probMin = probSP.min()
        probSP[probSP > 0.5 * probMin] = 0
        return probSP

    def matchScan(self, reading, estMovingDist, estMovingTheta, count, matchMax = True):
        """Iteratively find the best dx, dy and dtheta"""
        estimatedX, estimatedY, estimatedTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        rMeasure = np.asarray(rMeasure)
        if count == 1:
            return reading, 1
        # Coarse Search
        coarseSearchStep = self.coarseFactor * self.og.unitGridSize  # make this even number of unitGridSize for performance
        coarseSigma = self.scanSigmaInNumGrid / self.coarseFactor

        xRangeList, yRangeList, probSP = self.frameSearchSpace(estimatedX, estimatedY, coarseSearchStep, coarseSigma, self.missMatchProbAtCoarse)
        matchedPx, matchedPy, matchedReading, convTotal, coarseConfidence = self.searchToMatch(probSP, estimatedX, estimatedY,
            estimatedTheta, rMeasure, xRangeList, yRangeList, self.searchRadius,
                self.searchHalfRad, coarseSearchStep, estMovingDist, estMovingTheta,fineSearch=False, matchMax=matchMax)
        #########   For Debug Only  #############
        #if count > 14:
        #    self.plotMatchOverlay(probSP, matchedPx, matchedPy, matchedReading, xRangeList, yRangeList, coarseSearchStep)
        #########################################
        # Fine Search
        fineSearchStep = self.og.unitGridSize
        fineSigma = self.scanSigmaInNumGrid
        fineSearchHalfRad = self.searchHalfRad
        fineMissMatchProbAtFine = self.missMatchProbAtCoarse ** (2 / self.coarseFactor)
        xRangeList, yRangeList, probSP = self.frameSearchSpace(matchedReading['x'], matchedReading['y'], fineSearchStep, fineSigma, fineMissMatchProbAtFine)
        matchedPx, matchedPy, matchedReading, convTotal, fineConfidence = self.searchToMatch(probSP, matchedReading['x'],
            matchedReading['y'], matchedReading['theta'], matchedReading['range'], xRangeList, yRangeList,
                coarseSearchStep, fineSearchHalfRad, fineSearchStep, estMovingDist, estMovingTheta, fineSearch=True, matchMax=True)

        #########   For Debug Only  #############
        #if count > 0:
        #   self.plotMatchOverlay(probSP, matchedPx, matchedPy, matchedReading, xRangeList, yRangeList, fineSearchStep)
        #########################################
        return matchedReading, coarseConfidence

    def covertMeasureToXY(self, estimatedX, estimatedY, estimatedTheta, rMeasure):
        rads = np.linspace(estimatedTheta - self.og.lidarFOV / 2, estimatedTheta + self.og.lidarFOV / 2,
                           num=self.og.numSamplesPerRev)
        range_idx = rMeasure < self.og.lidarMaxRange
        rMeasureInRange = rMeasure[range_idx]
        rads = rads[range_idx]
        px = estimatedX + np.cos(rads) * rMeasureInRange
        py = estimatedY + np.sin(rads) * rMeasureInRange
        return px, py

    def searchToMatch(self, probSP, estimatedX, estimatedY, estimatedTheta, rMeasure, xRangeList, yRangeList,
                      searchRadius, searchHalfRad, unitLength, estMovingDist,  estMovingTheta, fineSearch = False, matchMax = True):
        px, py = self.covertMeasureToXY(estimatedX, estimatedY, estimatedTheta, rMeasure)
        numCellOfSearchRadius = int(searchRadius / unitLength)
        xMovingRange = np.arange(-numCellOfSearchRadius, numCellOfSearchRadius + 1)
        yMovingRange = np.arange(-numCellOfSearchRadius, numCellOfSearchRadius + 1)
        xv, yv = np.meshgrid(xMovingRange, yMovingRange)
        if fineSearch:
            rv, thetaWeight = np.zeros(xv.shape), np.zeros(xv.shape)
        else:
            rv = - (1 / (2 * self.moveRSigma ** 2)) * (np.sqrt((xv * unitLength) ** 2 + (yv * unitLength) ** 2) - (estMovingDist)) ** 2
            rrv = np.abs(np.sqrt((xv * unitLength) ** 2 + (yv * unitLength) ** 2) - estMovingDist)
            rv[rrv > self.maxMoveDeviation] = -100  # no points deviates more than maxMoveDeviation
            if estMovingTheta is not None:
                distv = np.sqrt(np.square(xv) + np.square(yv))
                distv[distv == 0] = 0.0001
                thetav = np.arccos((xv * math.cos(estMovingTheta) + yv * math.sin(estMovingTheta)) / distv)
                thetaWeight = -1 / (2 * self.turnSigma ** 2) * np.square(thetav)
            else:
                thetaWeight = np.zeros(xv.shape)

        xv = xv.reshape((xv.shape[0], xv.shape[1], 1))
        yv = yv.reshape((yv.shape[0], yv.shape[1], 1))
        thetaRange = np.arange(-searchHalfRad, searchHalfRad + self.og.angularStep, self.og.angularStep)
        convTotal = np.zeros((len(thetaRange), xv.shape[0], xv.shape[1]))
        for i, theta in enumerate(thetaRange):
            rotatedPx, rotatedPy = self.rotate((estimatedX, estimatedY), (px, py), theta)
            rotatedPxIdx, rotatedPyIdx = self.convertXYToSearchSpaceIdx(rotatedPx, rotatedPy, xRangeList[0],
                                                                        yRangeList[0], unitLength)
            uniqueRotatedPxPyIdx = np.unique(np.column_stack((rotatedPxIdx, rotatedPyIdx)), axis=0)
            rotatedPxIdx, rotatedPyIdx = uniqueRotatedPxPyIdx[:, 0], uniqueRotatedPxPyIdx[:, 1]
            #########   For Debug Only  #############
            # self.plotMatchOverlay(probSP, rotatedPx, rotatedPy, xRangeList, yRangeList, unitLength)
            #########################################
            rotatedPxIdx = rotatedPxIdx.reshape(1, 1, -1)
            rotatedPyIdx = rotatedPyIdx.reshape(1, 1, -1)
            rotatedPxIdx = rotatedPxIdx + xv
            rotatedPyIdx = rotatedPyIdx + yv
            convResult = probSP[rotatedPyIdx, rotatedPxIdx]
            convResultSum = np.sum(convResult, axis=2)
            convResultSum = convResultSum + rv + thetaWeight
            convTotal[i, :, :] = convResultSum
        if matchMax:
            maxIdx = np.unravel_index(convTotal.argmax(), convTotal.shape)
        else:
            convTotalFlatten = np.reshape(convTotal, -1)
            convTotalFlattenProb = np.exp(convTotalFlatten) / np.exp(convTotalFlatten).sum()
            maxIdx = np.random.choice(np.arange(convTotalFlatten.size), 1, p=convTotalFlattenProb)[0]
            maxIdx = np.unravel_index(maxIdx, convTotal.shape)

        confidence = np.sum(np.exp(convTotal))
        dx, dy, dtheta = xMovingRange[maxIdx[2]] * unitLength, yMovingRange[maxIdx[1]] * unitLength, thetaRange[maxIdx[0]]
        matchedReading = {"x": estimatedX + dx, "y": estimatedY + dy, "theta": estimatedTheta + dtheta,
                          "range": rMeasure}
        matchedPx, matchedPy = self.rotate((estimatedX, estimatedY), (px, py), dtheta)

        #########   For Debug Only  #############
        #self.plotMatchOverlay(probSP, matchedPx + dx, matchedPy + dy, matchedReading, xRangeList, yRangeList, unitLength)
        #########################################

        return matchedPx + dx, matchedPy + dy, matchedReading, convTotal, confidence

    def plotMatchOverlay(self, probSP, matchedPx, matchedPy, matchedReading, xRangeList, yRangeList, unitLength):
        plt.figure(figsize=(19.20, 19.20))
        plt.imshow(probSP, origin='lower')
        pxIdx, pyIdx = self.convertXYToSearchSpaceIdx(matchedPx, matchedPy, xRangeList[0], yRangeList[0], unitLength)
        plt.scatter(pxIdx, pyIdx, c='r', s=5)
        poseXIdx, poseYIdx = self.convertXYToSearchSpaceIdx(matchedReading['x'], matchedReading['y'], xRangeList[0], yRangeList[0], unitLength)
        plt.scatter(poseXIdx, poseYIdx, color='blue', s=50)
        plt.show()

    def rotate(self, origin, point, angle):
        """
        Rotate a point counterclockwise by a given angle around a given origin.
        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point
        qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
        qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
        return qx, qy

    def convertXYToSearchSpaceIdx(self, px, py, beginX, beginY, unitLength):
        xIdx = (((px - beginX) / unitLength)).astype(int)
        yIdx = (((py - beginY) / unitLength)).astype(int)
        return xIdx, yIdx

def updateEstimatedPose(currentRawReading, prevMatchedReading, prevRawReading, prevRawMovingTheta, prevMatchedMovingTheta):
    estimatedTheta = prevMatchedReading['theta'] + currentRawReading['theta'] - prevRawReading['theta']
    estimatedReading = {'x': prevMatchedReading['x'], 'y': prevMatchedReading['y'], 'theta': estimatedTheta, 'range': currentRawReading['range']}
    dx, dy = currentRawReading['x'] - prevRawReading['x'], currentRawReading['y'] - prevRawReading['y']
    estMovingDist = math.sqrt(dx**2 + dy**2)
    rawX, rawY, prevRawX, prevRawY = currentRawReading['x'], currentRawReading['y'], prevRawReading['x'], prevRawReading['y']
    rawXMove, rawYMove = rawX - prevRawX, rawY - prevRawY
    rawMove = math.sqrt((rawX - prevRawX) ** 2 + (rawY - prevRawY) ** 2)

    if rawMove > 0.3:
        if prevRawMovingTheta != None:
            if rawYMove > 0:
                rawMovingTheta = math.acos(rawXMove / rawMove)  # between -pi and +pi
            else:
                rawMovingTheta = -math.acos(rawXMove / rawMove)
            rawTurnTheta = rawMovingTheta - prevRawMovingTheta
            estMovingTheta = prevMatchedMovingTheta + rawTurnTheta
        else:
            if rawYMove > 0:
                rawMovingTheta = math.acos(rawXMove / rawMove)  # between -pi and +pi
            else:
                rawMovingTheta = -math.acos(rawXMove / rawMove)
            estMovingTheta = None
    else:
        rawMovingTheta = None
        estMovingTheta = None

    return estimatedReading, estMovingDist, estMovingTheta, rawMovingTheta

def updateTrajectory(matchedReading, xTrajectory, yTrajectory):
    x, y, theta, range = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
    xTrajectory.append(x)
    yTrajectory.append(y)

def getMovingTheta(matchedReading, xTrajectory, yTrajectory):
    x, y, theta, range = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
    prevX, prevY = xTrajectory[-1], yTrajectory[-1]
    xMove, yMove = x - prevX, y - prevY
    move = math.sqrt(xMove ** 2 + yMove ** 2)
    if move != 0:
        if yMove > 0:
            movingTheta = math.acos(xMove / move)
        else:
            movingTheta = -math.acos(xMove / move)
    else:
        movingTheta = None
    return movingTheta

def processSensorData(sensorData, og, sm):
    # gtData = readJson("../DataSet/PreprocessedData/intel_corrected_log") #########   For Debug Only  #############
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    colors = iter(cm.rainbow(np.linspace(1, 0, len(sensorData) + 1)))
    xTrajectory, yTrajectory = [], []
    for key in sorted(sensorData.keys()):
        count += 1
        print(count)
        if count == 1:
            #og.updateOccupancyGrid(sensorData[key])
            prevRawMovingTheta, prevMatchedMovingTheta = None, None
            matchedReading, confidence = sensorData[key], 1
            #prevGtReading = gtData[key]  #########   For Debug Only  #############
        else:
            currentRawReading = sensorData[key]
            estimatedReading, estMovingDist, estMovingTheta, rawMovingTheta = updateEstimatedPose(currentRawReading,
                prevMatchedReading, prevRawReading, prevRawMovingTheta, prevMatchedMovingTheta)
            matchedReading, confidence = sm.matchScan(estimatedReading, estMovingDist, estMovingTheta, count)
            prevRawMovingTheta = rawMovingTheta
            prevMatchedMovingTheta = getMovingTheta(matchedReading, xTrajectory, yTrajectory)
        #########   For Debug Only  #############
        #gtReading = gtData[key]
        #compareGT(currentRawReading, prevRawReading, matchedReading, prevMatchedReading, gtReading, prevGtReading)
        #prevGtReading = gtReading
        #########################################
        og.updateOccupancyGrid(matchedReading)
        updateTrajectory(matchedReading, xTrajectory, yTrajectory)
        prevMatchedReading, prevRawReading = matchedReading, sensorData[key]
        #if count == 100:
        #   break

    for i in range(len(xTrajectory)):
        plt.scatter(xTrajectory[i], yTrajectory[i], color=next(colors), s=35)
    plt.scatter(xTrajectory[0], yTrajectory[0], color='r', s=500)
    plt.scatter(xTrajectory[-1], yTrajectory[-1], color=next(colors), s=500)
    plt.plot(xTrajectory, yTrajectory)
    og.plotOccupancyGrid([-13, 20], [-25, 7], plotThreshold=False)

def readJson(jsonFile):
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        return input['map']

def compareGT(currentRawReading, prevRawReading, matchedReading, prevMatchedReading, gtReading, prevGtReading):
    print("true last pos x: " + str(prevGtReading['x']) + ", y: " + str(prevGtReading['y']))
    print("true curr pos x: " + str(gtReading['x']) + ", y: " + str(gtReading['y']))
    gtMoveX, gtMoveY = gtReading['x'] - prevGtReading['x'], gtReading['y'] - prevGtReading['y']
    print("true move x: " + str(gtMoveX) + ", y: " + str(gtMoveY) + ", r: " + str(
        math.sqrt(gtMoveX ** 2 + gtMoveY ** 2)))

    print("Estd last pos x: " + str(prevMatchedReading['x']) + ", y: " + str(prevMatchedReading['y']))
    print("Estd curr pos x: " + str(matchedReading['x']) + ", y: " + str(matchedReading['y']))
    rawEstMoveX, rawEstMoveY = currentRawReading['x'] - prevRawReading['x'], currentRawReading['y'] - \
                               prevRawReading['y']
    print("raw move x: " + str(rawEstMoveX) + ", y: " + str(rawEstMoveY) + ", r: " + str(
        math.sqrt(rawEstMoveX ** 2 + rawEstMoveY ** 2)))
    currMatchMinusRawMoveX = matchedReading['x'] - prevMatchedReading['x'] - rawEstMoveX
    currMatchMinusRawMoveY = matchedReading['y'] - prevMatchedReading['y'] - rawEstMoveY
    print(
        "compensate move x: " + str(currMatchMinusRawMoveX) + ", y: " + str(currMatchMinusRawMoveY) + ", r: " + str(
            math.sqrt(currMatchMinusRawMoveX ** 2 + currMatchMinusRawMoveY ** 2)))
    if math.sqrt(currMatchMinusRawMoveX ** 2 + currMatchMinusRawMoveY ** 2) > 1.5:
        a = 1

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, wallThickness, moveRSigma, maxMoveDeviation, turnSigma, \
        missMatchProbAtCoarse, coarseFactor = 1.4, 0.25, 2, 5 * unitGridSize, 0.1, 0.25, 0.3, 0.15, 5
    sensorData = readJson("../DataSet/PreprocessedData/intel_gfs")
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    initXY = sensorData[sorted(sensorData.keys())[0]]
    og = OccupancyGrid(initMapXLength, initMapYLength, initXY, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness)
    sm = ScanMatcher(og, scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, moveRSigma, maxMoveDeviation, turnSigma, missMatchProbAtCoarse, coarseFactor)
    processSensorData(sensorData, og, sm)

if __name__ == '__main__':
    main()