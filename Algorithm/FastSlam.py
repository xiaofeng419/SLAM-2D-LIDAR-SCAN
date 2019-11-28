import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from Utils.OccupancyGrid import OccupancyGrid
from Utils.ScanMatcher_OGBased import ScanMatcher
import math
import copy

class ParticleFilter:
    def __init__(self, numParticles, ogParameters, smParameters):
        self.numParticles = numParticles
        self.particles = []
        self.initParticles(ogParameters, smParameters)
        self.step = 0
        self.previousMatchedReading = None
        self.previousRawReading = None

    def initParticles(self, ogParameters, smParameters):
        for i in range(self.numParticles):
            p = Particle(ogParameters, smParameters)
            self.particles.append(p)

    def updateParticles(self, reading, count):
        for i in range(self.numParticles):
            self.particles[i].update(reading, count)

    def weightUnbalanced(self):
        self.normalizeWeights()
        variance = 0
        for i in range(self.numParticles):
            variance += (self.particles[i].weight - 1 / self.numParticles) ** 2
            #variance += self.particles[i].weight**2
        if variance > 0.88:
            return True
        else:
            return False

    def normalizeWeights(self):
        weightSum = 0
        for i in range(self.numParticles):
            weightSum += self.particles[i].weight
        for i in range(self.numParticles):
            self.particles[i].weight = self.particles[i].weight / weightSum

    def resample(self):

        # for particle in self.particles:
        #     particle.plotParticle()

        weights = np.zeros(self.numParticles)
        tempParticles = []
        for i in range(self.numParticles):
            weights[i] = self.particles[i].weight
            tempParticles.append(copy.deepcopy(self.particles[i]))
        resampledParticlesIdx = np.random.choice(np.arange(self.numParticles), self.numParticles, p=weights)
        for i in range(self.numParticles):
            self.particles[i] = copy.deepcopy(tempParticles[resampledParticlesIdx[i]])
            self.particles[i].weight = 1 / self.numParticles



class Particle:
    def __init__(self, ogParameters, smParameters):
        initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange, numSamplesPerRev = ogParameters
        scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, wallThickness, moveRSigma, missMatchProbAtCoarse, coarseFactor  = smParameters
        og = OccupancyGrid(initMapXLength, initMapYLength, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange, wallThickness)
        sm = ScanMatcher(og, scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, moveRSigma, missMatchProbAtCoarse, coarseFactor)
        self.og = og
        self.sm = sm
        self.xTrajectory = []
        self.yTrajectory = []
        self.weight = 1

    def updateEstimatedPose(self, currentRawReading):
        estimatedTheta = self.previousMatchedReading['theta'] + currentRawReading['theta'] - self.previousRawReading['theta']
        estimatedReading = {'x': self.previousMatchedReading['x'], 'y': self.previousMatchedReading['y'], 'theta': estimatedTheta,
                            'range': currentRawReading['range']}
        dx, dy = currentRawReading['x'] - self.previousRawReading['x'], currentRawReading['y'] - self.previousRawReading['y']
        estMovingDist = math.sqrt(dx ** 2 + dy ** 2)
        return estimatedReading, estMovingDist

    def update(self, reading, count):
        if count == 1:
            self.og.updateOccupancyGrid(reading)
            self.previousMatchedReading = reading
            self.previousRawReading = reading
        currentRawReading = reading
        estimatedReading, estMovingDist = self.updateEstimatedPose(currentRawReading)
        matchedReading, confidence = self.sm.matchScan(estimatedReading, estMovingDist, count, matchMax=False)
        self.updateTrajectoryPlot(matchedReading)
        self.og.updateOccupancyGrid(matchedReading)
        self.previousMatchedReading = matchedReading
        self.previousRawReading = reading
        self.weight *= confidence

    def updateTrajectoryPlot(self, matchedReading):
        x, y = matchedReading['x'], matchedReading['y']
        self.xTrajectory.append(x)
        self.yTrajectory.append(y)

    def plotParticle(self):
        plt.figure(figsize=(19.20, 19.20))
        plt.scatter(self.xTrajectory[0], self.yTrajectory[0], color='r', s=500)
        colors = iter(cm.rainbow(np.linspace(1, 0, len(self.xTrajectory) + 1)))
        for i in range(len(self.xTrajectory)):
            plt.scatter(self.xTrajectory[i], self.yTrajectory[i], color=next(colors), s=35)
        plt.scatter(self.xTrajectory[-1], self.yTrajectory[-1], color=next(colors), s=500)
        plt.plot(self.xTrajectory, self.yTrajectory)
        self.og.plotOccupancyGrid([-13, 20], [-25, 7], plotThreshold=False)

def processSensorData(pf, sensorData, plotTrajectory = True):
    # gtData = readJson("../DataSet/PreprocessedData/intel_corrected_log") #########   For Debug Only  #############
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    for key in sorted(sensorData.keys()):
        count += 1
        print(count)
        pf.updateParticles(sensorData[key], count)
        if pf.weightUnbalanced():
            pf.resample()
            print("resample")

        # if count == 100:
        #     break

    for particle in pf.particles:
        particle.plotParticle()

def readJson(jsonFile):
    with open(jsonFile, 'r') as f:
        input = json.load(f)
        return input['map']

def main():
    initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange = 10, 10, 0.02, np.pi, 10 # in Meters
    scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, wallThickness, moveRSigma, \
        missMatchProbAtCoarse, coarseFactor = 1.4, 0.25, 2, 5 * unitGridSize, 0.1, 0.2, 5
    sensorData = readJson("../DataSet/PreprocessedData/intel_gfs")
    numSamplesPerRev = len(sensorData[list(sensorData)[0]]['range'])  # Get how many points per revolution
    numParticles = 15
    ogParameters = [initMapXLength, initMapYLength, unitGridSize, lidarFOV, lidarMaxRange, numSamplesPerRev]
    smParameters = [scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, wallThickness, moveRSigma, missMatchProbAtCoarse, coarseFactor]
    pf = ParticleFilter(numParticles, ogParameters, smParameters)
    processSensorData(pf, sensorData, plotTrajectory=True)

if __name__ == '__main__':
    main()