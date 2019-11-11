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
        estimatedX, estimatedY, estimatedTheta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        rMeasure = np.asarray(rMeasure)
        xRangeList, yRangeList, extractedOG = self.extractLocalOG(estimatedX, estimatedY)
        probOG = gaussian_filter(extractedOG, sigma=self.scanSigmaInNumGrid)
        #plt.imshow(probOG)
        #plt.show()
        rads = np.linspace(estimatedTheta - self.og.lidarFOV / 2, estimatedTheta + self.og.lidarFOV / 2,
                           num=self.og.numSamplesPerRev)
        range_idx = rMeasure < self.og.lidarMaxRange
        rMeasure = rMeasure[range_idx]
        rads = rads[range_idx]
        px = estimatedX + np.cos(rads) * rMeasure
        py = estimatedY + np.sin(rads) * rMeasure
        xMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize).reshape(-1, 1)
        yMovingRange = np.arange(-self.searchRadius, self.searchRadius + self.og.unitGridSize, self.og.unitGridSize).reshape(-1, 1)
        xv, yv = np.meshgrid(xMovingRange, yMovingRange)
        xv = xv.reshape((xv.shape[0], xv.shape[1], 1))
        yv = yv.reshape((yv.shape[0], yv.shape[1], 1))
        maxMatch = 0
        for theta in np.arange(estimatedTheta - self.searchHalfRad, estimatedTheta + self.searchHalfRad + self.og.angularStep, self.og.angularStep):
            rotatedPx, rotatedPy = self.rotate((estimatedX, estimatedY), (px, py), theta)
            rotatedPx = rotatedPx.reshape(1, 1, -1)
            rotatedPy = rotatedPy.reshape(1, 1, -1)
            rotatedPx = rotatedPx + xv
            rotatedPy = rotatedPy + yv
            rotatedPxIdx, rotatedPyIdx = self.pxyToExtractedOGIdx(rotatedPx, rotatedPy, xRangeList, yRangeList)

            convResult = probOG[rotatedPxIdx, rotatedPyIdx]
            convResultMax = np.max(convResult, axis=2)
            if convResultMax.max() > maxMatch:
                maxMatch = convResultMax.max()
                maxIdx = np.unravel_index(convResultMax.argmax(), convResultMax.shape)
                maxTheta = theta
        dx = xMovingRange[maxIdx[1]]
        dy = yMovingRange[maxIdx[0]]
        matchedPx, matchedPy = self.rotate((estimatedX, estimatedY), (px, py), maxTheta)
        matchedPx = matchedPx + dx
        matchedPy = matchedPy + dy
        plt.scatter(matchedPx, matchedPy, c='r', s=35)
        self.og.plotOccupancyGrid()
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

        return xRangeList, yRangeList, self.og.occupancyGridVisited[yIdxList[0]: yIdxList[1], xIdxList[0]:xIdxList[1]] / \
            self.og.occupancyGridTotal[yIdxList[0]: yIdxList[1], xIdxList[0]:xIdxList[1]]


def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid = 0.5, 0.35, 2# 0.5m and 20deg
    wallThickness = 5 * unitGridSize
    jsonFile = "../DataSet/PreprocessedData/intel_corrected_log"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        sensorData = input['map']
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    spokesStartIdx = int(0) # theta= 0 is x direction. spokes=0 is y direction, the first ray of lidar scan direction. spokes increase clockwise
    og = OccupancyGrid(initMapXLength, initMapYLength, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness, spokesStartIdx)
    sm = ScanMatcher(og, scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid)
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    for key in sorted(sensorData.keys()):
        count += 1
        og.updateOccupancyGrid(sensorData[key])
        sm.matchScan(sensorData[key])
    map = np.flipud(1 - og.occupancyGridVisited / og.occupancyGridTotal)
    plt.matshow(map, cmap='gray', extent=[og.mapXLim[0], og.mapXLim[1], og.mapYLim[1], og.mapYLim[0]])
    plt.show()
    map = map < 0.5
    plt.matshow(map, cmap='gray', extent=[og.mapXLim[0], og.mapXLim[1], og.mapYLim[1], og.mapYLim[0]])
    plt.show()

if __name__ == '__main__':
    main()