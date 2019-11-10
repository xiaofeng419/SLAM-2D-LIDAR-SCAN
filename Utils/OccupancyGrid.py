import json
import numpy as np
import matplotlib.pyplot as plt

class occupancyGrid:
    def __init__(self, mapXLength, mapYLength, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness, spokesStartIdx):
        xNum = int(mapXLength / unitGridSize)
        yNum = int(mapYLength / unitGridSize)
        x = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=xNum + 1)
        y = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=yNum + 1)
        self.OccupancyGridX, self.OccupancyGridY = np.meshgrid(x, y)
        self.occupancyGridVisited = np.ones((xNum, yNum))
        self.occupancyGridTotal = np.ones((xNum, yNum))
        self.unitGridSize = unitGridSize
        self.lidarFOV = lidarFOV
        self.lidarMaxRange = lidarMaxRange
        self.wallThickness = wallThickness
        self.mapXLim = (self.OccupancyGridX[0, 0], self.OccupancyGridX[0, -1])
        self.mapYLim = (self.OccupancyGridY[0, 0], self.OccupancyGridY[-1, 0])
        self.numSamplesPerRev = numSamplesPerRev
        self.spokesStartIdx = spokesStartIdx
        self.angularStep = lidarFOV / numSamplesPerRev
        self.numSpokes = 2 * np.pi / self.angularStep
        xGrid, yGrid, bearingIdxGrid, rangeIdxGrid = self.spokesGrid()
        radByX, radByY, radByR = self.itemizeSpokesGrid(xGrid, yGrid, bearingIdxGrid, rangeIdxGrid)
        self.radByX = radByX
        self.radByY = radByY
        self.radByR = radByR

    def spokesGrid(self):
        """
        0th ray appear at north, then clock wise increase angle
        """
        numHalfElem = int(self.lidarMaxRange / self.unitGridSize)
        bearingIdxGrid = np.zeros((2 * numHalfElem + 1, 2 * numHalfElem + 1))
        x = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        y = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        xGrid, yGrid = np.meshgrid(x, y)
        bearingIdxGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] = np.rint((np.pi / 2 - np.arctan(
            yGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] / xGrid[:, numHalfElem + 1: 2 * numHalfElem + 1]))
                / np.pi / 2 * self.numSpokes - 0.5).astype(int)
        bearingIdxGrid[:, 0: numHalfElem] = np.fliplr(np.flipud(bearingIdxGrid))[:, 0: numHalfElem] + int(self.numSpokes / 2)
        bearingIdxGrid[0: numHalfElem, numHalfElem] = int(self.numSpokes / 2)
        rangeIdxGrid = np.sqrt(xGrid**2 + yGrid**2)
        return xGrid, yGrid, bearingIdxGrid, rangeIdxGrid

    def itemizeSpokesGrid(self, xGrid, yGrid, bearingIdxGrid, rangeIdxGrid):
        spokes = np.unique(bearingIdxGrid)
        radByX = []
        radByY = []
        radByR = []
        for i in spokes:
            idx = np.argwhere(bearingIdxGrid == i)
            radByX.append(xGrid[idx[:, 0], idx[:, 1]])
            radByY.append(yGrid[idx[:, 0], idx[:, 1]])
            radByR.append(rangeIdxGrid[idx[:, 0], idx[:, 1]])
        return radByX, radByY, radByR

    def expandOccupancyGridHelper(self, position, axis):
        gridShape = self.occupancyGridVisited.shape
        if axis == 0:
            insertion = np.ones((gridShape[0], int(gridShape[1] / 2)))
            if position == 0:
                x = self.OccupancyGridX[0]
                y = np.linspace(self.mapYLim[0] - int(gridShape[1] / 2) * self.unitGridSize, self.mapYLim[0],
                                num=int(gridShape[1] / 2), endpoint=False)
            else:
                x = self.OccupancyGridX[0]
                y = np.linspace(self.mapYLim[1] + self.unitGridSize, self.mapYLim[1] + (int(gridShape[1] / 2) + 1) * self.unitGridSize,
                                num=int(gridShape[1] / 2), endpoint=False)
        else:
            insertion = np.ones((int(gridShape[1]/2), gridShape[0]))
            if position == 0:
                y = self.OccupancyGridY[0]
                x = np.linspace(self.mapXLim[0] - int(gridShape[1] / 2) * self.unitGridSize, self.mapXLim[0],
                                num=int(gridShape[1] / 2), endpoint=False)
            else:
                y = self.OccupancyGridY[0]
                x = np.linspace(self.mapXLim[1] + self.unitGridSize, self.mapXLim[1] + (int(gridShape[1] / 2) + 1) * self.unitGridSize,
                                num=int(gridShape[1] / 2), endpoint=False)

        np.insert(self.occupancyGridVisited, position, insertion, axis=axis)
        np.insert(self.occupancyGridTotal, position, insertion, axis=axis)
        xv, yv = np.meshgrid(x, y)
        np.insert(self.OccupancyGridX, position, xv, axis=axis)
        np.insert(self.OccupancyGridY, position, yv, axis=axis)

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

    def convertRealXYToMapIdx(self, X, Y, mapXLim, mapYLim, unitGridSize):
        #mapXLim is (2,) array for left and right limit, same for mapYLim
        assert all(X > mapXLim[0] & X < mapXLim[1] & Y > mapYLim[0] & Y < mapYLim[1])
        xIdx = int(np.rint((X - mapXLim[0]) / unitGridSize + 1))
        yIdx = int(np.rint((Y - mapYLim[0]) / unitGridSize + 1))
        return xIdx, yIdx

    def checkMapToExpand(self, X, Y):
        if any(X < self.mapXLim[0]):
            return 1
        elif any(X > self.mapXLim[1]):
            return 2
        elif any(Y < self.mapYLim[0]):
            return 3
        elif any(Y > self.mapYLim[1]):
            return 4
        else:
            return -1

    def updateOccupancyGrid(self, reading):
        x, y, theta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        rMeasure = np.asarray(rMeasure)
        spokesOffsetIdxByTheta = int(np.rint((2 * np.pi - theta) / (2 * np.pi) * self.numSpokes))
        for i in range(self.numSamplesPerRev):
            idx = (self.spokesStartIdx + spokesOffsetIdxByTheta + i) % 360
            xInSpoke = self.radByX[idx]
            yInSpoke = self.radByY[idx]
            rInSpoke = self.radByR[idx]
            emptyIdx = np.argwhere(rInSpoke < rMeasure[i] - self.wallThickness / 2)
            occupiedIdx = np.argwhere((rInSpoke > rMeasure[i] - self.wallThickness / 2) & (rInSpoke < rMeasure[i] + self.wallThickness / 2))
            if len(emptyIdx) != 0 or len(occupiedIdx) != 0:
                expandDirection = self.checkMapToExpand(x + xInSpoke[emptyIdx], y + yInSpoke[emptyIdx])
                while (expandDirection != -1):
                    self.expandOccupancyGrid(expandDirection)
            #xIdx, yIdx = self.convertRealXYToMapIdx(x + xInSpoke[emptyIdx], y + yInSpoke[emptyIdx])

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 1, 1, 0.2, np.pi, 10 # in Meters
    wallThickness = 1.5 * unitGridSize
    jsonFile = "../DataSet/PreprocessedData/intel_corrected_log"
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        sensorData = input['map']
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    spokesStartIdx = int(0) # theta= 0 is x direction. spokes=0 is y direction, the first ray of lidar scan direction. spokes increase clockwise
    og = occupancyGrid(initMapXLength, initMapYLength, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness, spokesStartIdx)
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    for key in sorted(sensorData.keys()):
        count += 1
        og.updateOccupancyGrid(sensorData[key])

if __name__ == '__main__':
    main()