import json
import numpy as np
import matplotlib.pyplot as plt
from Utils.OccupancyGrid import OccupancyGrid
from scipy.ndimage import gaussian_filter

class ScanMatcher:
    def __init__(self, og, searchRadius, searchHalfRad, scanSigmaInNumGrid):
        self.searchRadius = searchRadius
        self.searchHalfRad = searchHalfRad
        self.og = og
        self.scanSigmaInNumGrid = scanSigmaInNumGrid
        self.lastScan = [] # 2D point clouds

    def generateProbSearchSpace(self, searchSpace):
        probSP = gaussian_filter(searchSpace, sigma=self.scanSigmaInNumGrid)
        maxCap = 1 / (2 * np.pi * self.scanSigmaInNumGrid ** 2)
        probSP[probSP > maxCap] = maxCap
        spShape = probSP.shape
        xx, yy = np.linspace(-1, 1, num=spShape[1]), np.linspace(-1, 1, num=spShape[0])
        xv, yv = np.meshgrid(xx, yy)
        rv = (np.square(xv) + np.square(yv))**1.5
        #probSP = probSP * rv
        return probSP

    def matchScan(self, reading):
        """Iteratively find the best dx, dy and dtheta"""
        estimatedX, estimatedY, estimatedTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        rMeasure = np.asarray(rMeasure)
        if self.lastScan == []:
            px, py = self.covertMeasureToXY(estimatedX, estimatedY, estimatedTheta, rMeasure)
            self.lastScan = [px, py]
            return reading
        xRangeList, yRangeList, searchSpace = self.frameSearchSpace(estimatedX, estimatedY)
        probSP = self.generateProbSearchSpace(searchSpace)
        matchedPx, matchedPy, matchedReading = self.searchToMatch(
            probSP, estimatedX, estimatedY, estimatedTheta, rMeasure, xRangeList, yRangeList)
        self.lastScan = [matchedPx, matchedPy]
        #########   For Debug Only  #############
        self.plotMatchOverlay(probSP, matchedPx, matchedPy, xRangeList, yRangeList)
        #########################################
        return matchedReading

    def covertMeasureToXY(self, estimatedX, estimatedY, estimatedTheta, rMeasure):
        rads = np.linspace(estimatedTheta - self.og.lidarFOV / 2, estimatedTheta + self.og.lidarFOV / 2,
                           num=self.og.numSamplesPerRev)
        range_idx = rMeasure < self.og.lidarMaxRange
        rMeasureInRange = rMeasure[range_idx]
        rads = rads[range_idx]
        px = estimatedX + np.cos(rads) * rMeasureInRange
        py = estimatedY + np.sin(rads) * rMeasureInRange
        return px, py

    def searchToMatch(self, probSP, estimatedX, estimatedY, estimatedTheta, rMeasure, xRangeList, yRangeList):
        px, py = self.covertMeasureToXY(estimatedX, estimatedY, estimatedTheta, rMeasure)
        xMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize)
        yMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize)
        xv, yv = np.meshgrid(xMovingRange, yMovingRange)
        xv = xv.reshape((xv.shape[0], xv.shape[1], 1))
        yv = yv.reshape((yv.shape[0], yv.shape[1], 1))
        maxMatchScore, maxIdx = 0, None
        for theta in np.arange(-self.searchHalfRad, self.searchHalfRad + self.og.angularStep, self.og.angularStep):
            rotatedPx, rotatedPy = self.rotate((estimatedX, estimatedY), (px, py), theta)
            #########   For Debug Only  #############
            #self.plotMatchOverlay(probSP, rotatedPx, rotatedPy, xRangeList, yRangeList)
            #########################################
            rotatedPx = rotatedPx.reshape(1, 1, -1)
            rotatedPy = rotatedPy.reshape(1, 1, -1)
            rotatedPx = rotatedPx + xv
            rotatedPy = rotatedPy + yv
            rotatedPxIdx, rotatedPyIdx = self.convertXYToSearchSpaceIdx(rotatedPx, rotatedPy, xRangeList[0], yRangeList[0])
            convResult = probSP[rotatedPyIdx, rotatedPxIdx]
            convResultSum = np.sum(convResult, axis=2)
            if convResultSum.max() > maxMatchScore:
                maxMatchScore = convResultSum.max()
                maxIdx = np.unravel_index(convResultSum.argmax(), convResultSum.shape)
                dTheta = theta
                #########   For Debug Only  #############
                # matchedPx, matchedPy = self.rotate((estimatedX, estimatedY), (px, py), dTheta)
                # dx, dy = xMovingRange[maxIdx[1]], yMovingRange[maxIdx[0]]
                # self.plotMatchOverlay(probSP, matchedPx + dx, matchedPy + dy, xRangeList, yRangeList)
                # a = 1
                #########################################
        if maxIdx is None:
            dx, dy, dtheta = 0, 0, 0
        else:
            dx, dy = xMovingRange[maxIdx[1]], yMovingRange[maxIdx[0]]
        matchedReading = {"x": estimatedX + dx, "y": estimatedY + dy, "theta": estimatedTheta + dTheta,
                          "range": rMeasure}
        matchedPx, matchedPy = self.rotate((estimatedX, estimatedY), (px, py), dTheta)
        return matchedPx + dx, matchedPy + dy, matchedReading

    def plotMatchOverlay(self, probSP, matchedPx, matchedPy, xRangeList, yRangeList):
        plt.figure(figsize=(19.20, 19.20))
        plt.imshow(probSP)
        pxIdx, pyIdx = self.convertXYToSearchSpaceIdx(matchedPx, matchedPy, xRangeList[0], yRangeList[0])
        plt.scatter(pxIdx, pyIdx, c='r', s=1)
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

    def convertXYToSearchSpaceIdx(self, px, py, beginX, beginY):
        xIdx = (((px - beginX) / self.og.unitGridSize)).astype(int)
        yIdx = (((py - beginY) / self.og.unitGridSize)).astype(int)
        return xIdx, yIdx

    def frameSearchSpace(self, estimatedX, estimatedY):
        maxScanRadius = 1.1 * self.og.lidarMaxRange + self.searchRadius
        xRangeList = [estimatedX - maxScanRadius, estimatedX + maxScanRadius]
        yRangeList = [estimatedY - maxScanRadius, estimatedY + maxScanRadius]
        idxEndX, idxEndY = int((xRangeList[1] - xRangeList[0]) / self.og.unitGridSize),  int((yRangeList[1] - yRangeList[0]) / self.og.unitGridSize)
        searchSpace = np.zeros((idxEndY + 1, idxEndX + 1))
        lastScanIdx = self.convertXYToSearchSpaceIdx(self.lastScan[0], self.lastScan[1], xRangeList[0], yRangeList[0])
        searchSpace[lastScanIdx[1], lastScanIdx[0]] = 1
        return xRangeList, yRangeList, searchSpace

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid = 2.2, 0.35, 3# 0.5m and 20deg
    wallThickness = 5 * unitGridSize
    jsonFile = "../DataSet/PreprocessedData/intel_gfs"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        sensorData = input['map']
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    spokesStartIdx = int(0) # theta= 0 is x direction. spokes=0 is -y direction, the first ray of lidar scan direction. spokes increase counter-clockwise
    og = OccupancyGrid(initMapXLength, initMapYLength, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness, spokesStartIdx)
    sm = ScanMatcher(og, scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid)
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    for key in sorted(sensorData.keys()):
        count += 1
        if count == 1:
            og.updateOccupancyGrid(sensorData[key])
            previousMatchedReading = sensorData[key]
            previousRawReading = sensorData[key]

        currentRawReading = sensorData[key]
        estimatedX = previousMatchedReading['x'] + currentRawReading['x'] - previousRawReading['x']
        estimatedY = previousMatchedReading['y'] + currentRawReading['y'] - previousRawReading['y']
        estimatedTheta = previousMatchedReading['theta'] + currentRawReading['theta'] - previousRawReading['theta']
        estimatedReading = {'x': estimatedX, 'y': estimatedY, 'theta': estimatedTheta, 'range': sensorData[key]['range']}
        matchedReading = sm.matchScan(estimatedReading)
        og.updateOccupancyGrid(matchedReading)
        #og.plotOccupancyGrid(plotThreshold=False)
        previousMatchedReading = matchedReading
        previousRawReading = sensorData[key]
        print(count)
        if count == 100:
            break
    og.plotOccupancyGrid(plotThreshold=False)

if __name__ == '__main__':
    main()