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

    def matchScan(self, reading):
        """Iteratively find the best dx, dy and dtheta"""
        estimatedX, estimatedY, estimatedTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        rMeasure = np.asarray(rMeasure)
        xRangeList, yRangeList, extractedOG = self.extractLocalOG(estimatedX, estimatedY)
        probOG = gaussian_filter(extractedOG, sigma=self.scanSigmaInNumGrid)
        maxCap = 1 / np.sqrt(2 * np.pi * self.scanSigmaInNumGrid**2)
        probOG[probOG > maxCap] = maxCap
        rads = np.linspace(estimatedTheta - self.og.lidarFOV / 2, estimatedTheta + self.og.lidarFOV / 2,
                           num=self.og.numSamplesPerRev)
        range_idx = rMeasure < self.og.lidarMaxRange
        rMeasureInRange = rMeasure[range_idx]
        rads = rads[range_idx]
        px = estimatedX + np.cos(rads) * rMeasureInRange
        py = estimatedY + np.sin(rads) * rMeasureInRange
        xMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize)
        yMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize)
        xv, yv = np.meshgrid(xMovingRange, yMovingRange)
        xv = xv.reshape((xv.shape[0], xv.shape[1], 1))
        yv = yv.reshape((yv.shape[0], yv.shape[1], 1))
        maxMatch = 0
        for theta in np.arange(-self.searchHalfRad, self.searchHalfRad + self.og.angularStep, self.og.angularStep):
            rotatedPx, rotatedPy = self.rotate((estimatedX, estimatedY), (px, py), theta)
            rotatedPx = rotatedPx.reshape(1, 1, -1)
            rotatedPy = rotatedPy.reshape(1, 1, -1)
            rotatedPx = rotatedPx + xv
            rotatedPy = rotatedPy + yv
            rotatedPxIdx, rotatedPyIdx = self.pxyToExtractedOGIdx(rotatedPx, rotatedPy, xRangeList, yRangeList)
            convResult = probOG[rotatedPyIdx, rotatedPxIdx]
            convResultSum = np.sum(convResult, axis=2)
            if convResultSum.max() > maxMatch:
                maxMatch = convResultSum.max()
                maxIdx = np.unravel_index(convResultSum.argmax(), convResultSum.shape)
                dTheta = theta
        dx = xMovingRange[maxIdx[1]]
        dy = yMovingRange[maxIdx[0]]
        matchedReading = {"x": estimatedX + dx, "y": estimatedY + dy, "theta": estimatedTheta + dTheta, "range": rMeasure}

        plt.figure(figsize=(19.20, 19.20))
        matchedPx, matchedPy = self.rotate((estimatedX, estimatedY), (px, py), dTheta)
        matchedPx = matchedPx + dx
        matchedPy = matchedPy + dy
        plt.imshow(extractedOG)
        pxIdx, pyIdx = self.pxyToExtractedOGIdx(matchedPx, matchedPy, xRangeList, yRangeList)
        plt.scatter(pxIdx, pyIdx, c='r', s=1)
        plt.show()
        return matchedReading

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

    def pxyToExtractedOGIdx(self, px, py, xRangeList, yRangeList):
        xIdx = (np.rint((px - xRangeList[0]) / self.og.unitGridSize)).astype(int)
        yIdx = (np.rint((py - yRangeList[0]) / self.og.unitGridSize)).astype(int)
        return xIdx, yIdx

    def extractLocalOG(self, estimatedX, estimatedY):
        maxScanRadius = self.og.lidarMaxRange + self.searchRadius
        xRangeList = [estimatedX - maxScanRadius, estimatedX + maxScanRadius]
        yRangeList = [estimatedY - maxScanRadius, estimatedY + maxScanRadius]
        self.og.checkAndExapndOG(xRangeList, yRangeList)
        xIdxList, yIdxList = self.og.convertRealXYToMapIdx(xRangeList, yRangeList)
        extractedOg = self.og.occupancyGridVisited[yIdxList[0]: yIdxList[1], xIdxList[0]:xIdxList[1]] / \
            self.og.occupancyGridTotal[yIdxList[0]: yIdxList[1], xIdxList[0]:xIdxList[1]]
        extractedOg[extractedOg > 0.5] = 1
        extractedOg[extractedOg <= 0.5] = 0
        return xRangeList, yRangeList, extractedOg

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid = 0.5, 0.35, 5# 0.5m and 20deg
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
    og.plotOccupancyGrid(plotThreshold=False)

if __name__ == '__main__':
    main()