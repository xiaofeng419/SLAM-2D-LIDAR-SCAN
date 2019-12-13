import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class OccupancyGrid:
    def __init__(self, mapXLength, mapYLength, initXY, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness):
        xNum = int(mapXLength / unitGridSize)
        yNum = int(mapYLength / unitGridSize)
        x = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=xNum + 1) + initXY['x']
        y = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=yNum + 1) + initXY['y']
        self.OccupancyGridX, self.OccupancyGridY = np.meshgrid(x, y)
        self.occupancyGridVisited = np.ones((xNum + 1, yNum + 1))
        self.occupancyGridTotal = 2 * np.ones((xNum + 1, yNum + 1))
        self.unitGridSize = unitGridSize
        self.lidarFOV = lidarFOV
        self.lidarMaxRange = lidarMaxRange
        self.wallThickness = wallThickness
        self.mapXLim = [self.OccupancyGridX[0, 0], self.OccupancyGridX[0, -1]]
        self.mapYLim = [self.OccupancyGridY[0, 0], self.OccupancyGridY[-1, 0]]
        self.numSamplesPerRev = numSamplesPerRev
        self.angularStep = lidarFOV / numSamplesPerRev
        self.numSpokes = int(np.rint(2 * np.pi / self.angularStep))
        xGrid, yGrid, bearingIdxGrid, rangeIdxGrid = self.spokesGrid()
        radByX, radByY, radByR = self.itemizeSpokesGrid(xGrid, yGrid, bearingIdxGrid, rangeIdxGrid)
        self.radByX = radByX
        self.radByY = radByY
        self.radByR = radByR
        # theta= 0 is x direction. spokes=0 is y direc tion, spokesStartIdx is the first ray of lidar scan direction. spokes increase counter-clockwise
        self.spokesStartIdx = int(((self.numSpokes / 2 - self.numSamplesPerRev) / 2) % self.numSpokes)

    def spokesGrid(self):
        # 0th ray is at south, then counter-clock wise increases. Theta 0 is at east.
        numHalfElem = int(self.lidarMaxRange / self.unitGridSize)
        bearingIdxGrid = np.zeros((2 * numHalfElem + 1, 2 * numHalfElem + 1))
        x = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        y = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        xGrid, yGrid = np.meshgrid(x, y)
        bearingIdxGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] = np.rint((np.pi / 2 + np.arctan(
            yGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] / xGrid[:, numHalfElem + 1: 2 * numHalfElem + 1]))
                / np.pi / 2 * self.numSpokes - 0.5).astype(int)
        bearingIdxGrid[:, 0: numHalfElem] = np.fliplr(np.flipud(bearingIdxGrid))[:, 0: numHalfElem] + int(self.numSpokes / 2)
        bearingIdxGrid[numHalfElem + 1: 2 * numHalfElem + 1, numHalfElem] = int(self.numSpokes / 2)
        rangeIdxGrid = np.sqrt(xGrid**2 + yGrid**2)
        return xGrid, yGrid, bearingIdxGrid, rangeIdxGrid

    def itemizeSpokesGrid(self, xGrid, yGrid, bearingIdxGrid, rangeIdxGrid):
        # Due to discretization, later theta added could lead to up to 1 deg discretization error
        radByX = []
        radByY = []
        radByR = []
        for i in range(self.numSpokes):
            idx = np.argwhere(bearingIdxGrid == i)
            radByX.append(xGrid[idx[:, 0], idx[:, 1]])
            radByY.append(yGrid[idx[:, 0], idx[:, 1]])
            radByR.append(rangeIdxGrid[idx[:, 0], idx[:, 1]])
        return radByX, radByY, radByR

    def expandOccupancyGridHelper(self, position, axis):
        gridShape = self.occupancyGridVisited.shape
        if axis == 0:
            insertion = np.ones((int(gridShape[0] / 5),  gridShape[1]))
            if position == 0:
                x = self.OccupancyGridX[0]
                y = np.linspace(self.mapYLim[0] - int(gridShape[0] / 5) * self.unitGridSize, self.mapYLim[0],
                                num=int(gridShape[0] / 5), endpoint=False)
            else:
                x = self.OccupancyGridX[0]
                y = np.linspace(self.mapYLim[1] + self.unitGridSize, self.mapYLim[1] + (int(gridShape[0] / 5) ) * self.unitGridSize,
                                num=int(gridShape[0] / 5), endpoint=False)
        else:
            insertion = np.ones((gridShape[0], int(gridShape[1] / 5)))
            if position == 0:
                y = self.OccupancyGridY[:, 0]
                x = np.linspace(self.mapXLim[0] - int(gridShape[1] / 5) * self.unitGridSize, self.mapXLim[0],
                                num=int(gridShape[1] / 5), endpoint=False)
            else:
                y = self.OccupancyGridY[:, 0]
                x = np.linspace(self.mapXLim[1] + self.unitGridSize, self.mapXLim[1] + (int(gridShape[1] / 5)) * self.unitGridSize,
                                num=int(gridShape[1] / 5), endpoint=False)
        self.occupancyGridVisited = np.insert(self.occupancyGridVisited, [position], insertion, axis=axis)
        self.occupancyGridTotal = np.insert(self.occupancyGridTotal, [position], 2 * insertion, axis=axis)
        xv, yv = np.meshgrid(x, y)
        self.OccupancyGridX = np.insert(self.OccupancyGridX, [position], xv, axis=axis)
        self.OccupancyGridY = np.insert(self.OccupancyGridY, [position], yv, axis=axis)
        self.mapXLim[0] = self.OccupancyGridX[0, 0]
        self.mapXLim[1] = self.OccupancyGridX[0, -1]
        self.mapYLim[0] = self.OccupancyGridY[0, 0]
        self.mapYLim[1] = self.OccupancyGridY[-1, 0]

    def expandOccupancyGrid(self, expandDirection):
        gridShape = self.occupancyGridVisited.shape
        if expandDirection == 1:
            self.expandOccupancyGridHelper(0, 1)
        elif expandDirection == 2:
            self.expandOccupancyGridHelper(gridShape[1], 1)
        elif expandDirection == 3:
            self.expandOccupancyGridHelper(0, 0)
        else:
            self.expandOccupancyGridHelper(gridShape[0], 0)

    def convertRealXYToMapIdx(self, x, y):
        #mapXLim is (2,) array for left and right limit, same for mapYLim
        xIdx = (np.rint((x - self.mapXLim[0]) / self.unitGridSize)).astype(int)
        yIdx = (np.rint((y - self.mapYLim[0]) / self.unitGridSize)).astype(int)
        return xIdx, yIdx

    def checkMapToExpand(self, x, y):
        if any(x < self.mapXLim[0]):
            return 1
        elif any(x > self.mapXLim[1]):
            return 2
        elif any(y < self.mapYLim[0]):
            return 3
        elif any(y > self.mapYLim[1]):
            return 4
        else:
            return -1

    def checkAndExapndOG(self, x, y):
        """check x, y (vector points) are inside OG. If not, expand OG."""
        expandDirection = self.checkMapToExpand(x, y)
        while (expandDirection != -1):
            self.expandOccupancyGrid(expandDirection)
            expandDirection = self.checkMapToExpand(x, y)

    def updateOccupancyGrid(self, reading, dTheta = 0, update=True):
        x, y, theta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        theta += dTheta
        rMeasure = np.asarray(rMeasure)
        spokesOffsetIdxByTheta = int(np.rint(theta / (2 * np.pi) * self.numSpokes))
        emptyXList, emptyYList, occupiedXList, occupiedYList = [], [], [], []
        for i in range(self.numSamplesPerRev):
            spokeIdx = int(np.rint((self.spokesStartIdx + spokesOffsetIdxByTheta + i) % self.numSpokes))
            xAtSpokeDir = self.radByX[spokeIdx]
            yAtSpokeDir = self.radByY[spokeIdx]
            rAtSpokeDir = self.radByR[spokeIdx]
            if rMeasure[i] < self.lidarMaxRange:
                emptyIdx = np.argwhere(rAtSpokeDir < rMeasure[i] - self.wallThickness / 2)
            else:
                emptyIdx = []
            occupiedIdx = np.argwhere(
                (rAtSpokeDir > rMeasure[i] - self.wallThickness / 2) & (rAtSpokeDir < rMeasure[i] + self.wallThickness / 2))
            xEmptyIdx, yEmptyIdx = self.convertRealXYToMapIdx(x + xAtSpokeDir[emptyIdx], y + yAtSpokeDir[emptyIdx])
            xOccupiedIdx, yOccupiedIdx = self.convertRealXYToMapIdx(x + xAtSpokeDir[occupiedIdx], y + yAtSpokeDir[occupiedIdx])
            if update:
                self.checkAndExapndOG(x + xAtSpokeDir[occupiedIdx], y + yAtSpokeDir[occupiedIdx])
                if len(emptyIdx) != 0:
                    self.occupancyGridTotal[yEmptyIdx, xEmptyIdx] += 1
                if len(occupiedIdx) != 0:
                    self.occupancyGridVisited[yOccupiedIdx, xOccupiedIdx] += 2
                    self.occupancyGridTotal[yOccupiedIdx, xOccupiedIdx] += 2
            else:
                emptyXList.extend(x + xAtSpokeDir[emptyIdx])
                emptyYList.extend(y + yAtSpokeDir[emptyIdx])
                occupiedXList.extend(x + xAtSpokeDir[occupiedIdx])
                occupiedYList.extend(y + yAtSpokeDir[occupiedIdx])
        if not update:
            return np.asarray(emptyXList), np.asarray(emptyYList), np.asarray(occupiedXList), np.asarray(occupiedYList)

    def plotOccupancyGrid(self, xRange = None, yRange= None, plotThreshold = True):
        if xRange is None or xRange[0] < self.mapXLim[0] or xRange[1] > self.mapXLim[1]:
            xRange = self.mapXLim
        if yRange is None or yRange[0] < self.mapYLim[0] or yRange[1] > self.mapYLim[1]:
            yRange = self.mapYLim
        ogMap = self.occupancyGridVisited / self.occupancyGridTotal
        xIdx, yIdx = self.convertRealXYToMapIdx(xRange, yRange)
        ogMap = ogMap[yIdx[0]: yIdx[1], xIdx[0]: xIdx[1]]
        ogMap = np.flipud(1 - ogMap)
        plt.imshow(ogMap, cmap='gray', extent=[xRange[0], xRange[1], yRange[0], yRange[1]])
        plt.show()
        if plotThreshold:
            ogMap = ogMap >= 0.5
            plt.matshow(ogMap, cmap='gray', extent=[xRange[0], xRange[1], yRange[0], yRange[1]])
            plt.show()

def updateTrajectoryPlot(matchedReading, xTrajectory, yTrajectory, colors, count):
    x, y, theta, range = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
    xTrajectory.append(x)
    yTrajectory.append(y)
    if count % 1 == 0:
        plt.scatter(x, y, color=next(colors), s=35)

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    wallThickness = 7 * unitGridSize
    jsonFile = "../DataSet/PreprocessedData/intel_gfs"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        sensorData = input['map']
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    initXY = sensorData[sorted(sensorData.keys())[0]]
    og = OccupancyGrid(initMapXLength, initMapYLength, initXY, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness)
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    xTrajectory, yTrajectory = [], []
    colors = iter(cm.rainbow(np.linspace(1, 0, len(sensorData) + 1)))
    for key in sorted(sensorData.keys()):
        count += 1
        og.updateOccupancyGrid(sensorData[key])
        updateTrajectoryPlot(sensorData[key], xTrajectory, yTrajectory, colors, count)
        #if count == 100:
        #   break

    plt.scatter(xTrajectory[0], yTrajectory[0], color='r', s=500)
    plt.scatter(xTrajectory[-1], yTrajectory[-1], color=next(colors), s=500)
    plt.plot(xTrajectory, yTrajectory)
    #og.plotOccupancyGrid([-12, 20], [-23.5, 7])
    #og.plotOccupancyGrid(xRange =[-11.9, 20], yRange =[-23.5, 6.    og.plotOccupancyGrid()

if __name__ == '__main__':
    main()