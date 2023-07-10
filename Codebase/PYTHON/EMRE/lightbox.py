###########################################################################################
# Lightbox Simulation Code
# Rewritten by Emre Mengi from original code by Tarek Zohdi.
# Copyright 2023 Tarek Zohdi, Emre Mengi. All rights reserved. 
###########################################################################################

#%% Importing Packages

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from mpl_toolkits.mplot3d import Axes3D
from copy import deepcopy
from matplotlib import animation
from matplotlib import rc
from numpy.random import uniform
import os
import csv 

# %% Simulation Constants

# Simulation constants
lightboxConstants = {

    # time-stepping constants
    'stepFactor': 2.0, # factor by which to increase time step size
    'timeMultiplier': 5.0, # multiplier to determine final time

    # raytracing constants
    'c': 3e8, # speed of light (m/s)
    'rayDensity': 200, # number of rays per face
    'reflections': True, # marks if reflections are allowed
    'rayTol': 0.000001, # threshold for ray absorption

    'X1MinusOn': True, # marks if -x1 face will have rays
    'X1PlusOn': True, # marks if +x1 face will have rays
    'X2MinusOn': True, # marks if -x2 face will have rays
    'X2PlusOn': True, # marks if +x2 face will have rays
    'X3MinusOn': True, # marks if -x3 face will have rays
    'X3PlusOn': True, # marks if +x3 face will have rays

    # voxel constants
    'deltaX1': 0.05, # voxel length (m)
    'deltaX2': 0.05, # voxel width (m)
    'deltaX3': 0.05, # voxel height (m)
    'nVoxelX1': 20, # number of voxels in x1 direction
    'nVoxelX2': 20, # number of voxels in x2 direction
    'nVoxelX3': 20, # number of voxels in x3 direction

    # target constants
    # 'nTargetsX1': 9, # number of targets in x1 direction
    # 'nTargetsX2': 3, # number of targets in x2 direction
    # 'nTargetsX3': 3, # number of targets in x3 direction
    # 'rackX1X3': True, # marks if rack is in x1-x3 plane
    # 'rackX2X3': True, # marks if rack is in x2-x3 plane
    # 'rackX1X2': True, # marks if rack is in x1-x2 plane
    # 'rackX1': False, # marks if rack is in x1 direction
    # 'rackX2': False, # marks if rack is in x2 direction
    # 'rackX3': False, # marks if rack is in x3 direction

    'nTargetsX1': 9, # number of targets in x1 direction
    'nTargetsX2': 3, # number of targets in x2 direction
    'nTargetsX3': 3, # number of targets in x3 direction
    'rackX1X3': False, # marks if rack is in x1-x3 plane
    'rackX2X3': False, # marks if rack is in x2-x3 plane
    'rackX1X2': False, # marks if rack is in x1-x2 plane
    'rackX1': True, # marks if rack is in x1 direction
    'rackX2': False, # marks if rack is in x2 direction
    'rackX3': False, # marks if rack is in x3 direction

    # power constants
    'powerMagnitude': 1.0, # magnitude of power (to be multiplied by power dependency parameters)

    # misc constants
    'randseed': np.random.RandomState(201) # seed for ray variable randomizer

}

# Movie settings
movieSettings = {
    
    # toggle what to save/plot
    'targetSave': True, # save target data
    'raySave': True, # save ray data

    'movieFrames': 100, # number of frames in movie
    'movieTime': 0.0 # start time of movie [s]

}
    
#%% Defining Functions

# Lightbox simulation function
def lightboxSim(Lam, consts, movieSettings):

    # print line to indicate start of simulation
    print('----------------------------------------')
    print('Simulation started.')

    ########## Assigning simulation constants to associated variables ##########
    # time-stepping constants
    stepFactor = consts['stepFactor'] # factor by which to increase time step size
    timeMultiplier = consts['timeMultiplier'] # multiplier to determine final time
    # raytracing constants
    c = consts['c'] # speed of light (m/s)
    rayDensity = consts['rayDensity'] # number of rays per face
    reflections = consts['reflections'] # marks if reflections are allowed
    X1MinusOn = consts['X1MinusOn'] # marks if -x1 face will have rays
    rayTol = consts['rayTol'] # threshold for ray absorption

    X1PlusOn = consts['X1PlusOn'] # marks if +x1 face will have rays
    X2MinusOn = consts['X2MinusOn'] # marks if -x2 face will have rays
    X2PlusOn = consts['X2PlusOn'] # marks if +x2 face will have rays
    X3MinusOn = consts['X3MinusOn'] # marks if -x3 face will have rays
    X3PlusOn = consts['X3PlusOn'] # marks if +x3 face will have rays
    # voxel constants
    deltaX1 = consts['deltaX1'] # voxel length (m)
    deltaX2 = consts['deltaX2'] # voxel width (m)
    deltaX3 = consts['deltaX3'] # voxel height (m)
    nVoxelX1 = consts['nVoxelX1'] # number of voxels in x1 direction
    nVoxelX2 = consts['nVoxelX2'] # number of voxels in x2 direction
    nVoxelX3 = consts['nVoxelX3'] # number of voxels in x3 direction
    # target constants
    nTargetsX1 = consts['nTargetsX1'] # number of targets in x1 direction
    nTargetsX2 = consts['nTargetsX2'] # number of targets in x2 direction
    nTargetsX3 = consts['nTargetsX3'] # number of targets in x3 direction
    rackX1X3 = consts['rackX1X3'] # marks if rack is in x1-x3 plane
    rackX2X3 = consts['rackX2X3'] # marks if rack is in x2-x3 plane
    rackX1X2 = consts['rackX1X2'] # marks if rack is in x1-x2 plane
    rackX1 = consts['rackX1'] # marks if rack is in x1 direction
    rackX2 = consts['rackX2'] # marks if rack is in x2 direction
    rackX3 = consts['rackX3'] # marks if rack is in x3 direction
    # power constants
    powerMagnitude = consts['powerMagnitude'] # magnitude of power (to be multiplied by power dependency parameters)
    # misc constants
    randseed = consts['randseed'] # seed for ray variable randomizer

    ########### Assigning design string to associated variables ###############  
    # beam spread parameters (18 total)
    A1 = Lam[0]
    A2 = Lam[1]
    A3 = Lam[2]
    A4 = Lam[3]
    A5 = Lam[4]
    A6 = Lam[5]
    A7 = Lam[6]
    A8 = Lam[7]
    A9 = Lam[8]
    A10 = Lam[9]
    A11 = Lam[10]
    A12 = Lam[11]
    A13 = Lam[12]
    A14 = Lam[13]
    A15 = Lam[14]
    A16 = Lam[15]
    A17 = Lam[16]
    A18 = Lam[17]
    # sourcetube parameters (12 total)
    forPlusX3SourcetubeX1 = Lam[18]
    forPlusX3SourcetubeX2 = Lam[19]
    forMinusX3SourcetubeX1 = Lam[20]
    forMinusX3SourcetubeX2 = Lam[21]

    forPlusX2SourcetubeX1 = Lam[22]
    forPlusX2SourcetubeX3 = Lam[23]
    forMinusX2SourcetubeX1 = Lam[24]
    forMinusX2SourcetubeX3 = Lam[25]

    forPlusX1SourcetubeX2 = Lam[26]
    forPlusX1SourcetubeX3 = Lam[27]
    forMinusX1SourcetubeX2 = Lam[28]
    forMinusX1SourcetubeX3 = Lam[29]
    # power dependency parameters (6 total)
    magFacePowerX1Minus = Lam[30]
    magFacePowerX1Plus = Lam[31]
    magFacePowerX2Minus = Lam[32]
    magFacePowerX2Plus = Lam[33]
    magFacePowerX3Minus = Lam[34]
    magFacePowerX3Plus = Lam[35]

    ############################## Extracting movie settings ##################

    # toggle what to save/plot
    targetSave = movieSettings['targetSave'] # save target data
    raySave = movieSettings['raySave'] # save ray data

    movieFrames = movieSettings['movieFrames'] # number of frames in movie
    movieTime = movieSettings['movieTime'] # start time of movie [s]

    if raySave == True: # if ray data is to be saved  
            # delete data file from previous run
            if os.path.exists("ray.dat"):
                os.remove("ray.dat")
                # Print the statement once
                # the file is deleted   
                print("Previous file deleted.")                         
    if targetSave == True: # if target data is to be saved
        # delete data file from previous run
        if os.path.exists("target.dat"):
            os.remove("target.dat")
            # Print the statement once
            # the file is deleted 
            print("Previous file deleted.") 

    #################### Pre-allocation of variables ##########################

    ## target variables
    # target locations
    targetX1 = np.empty(1000)
    targetX2 = np.empty(1000)
    targetX3 = np.empty(1000)
    # target radii
    targetRadiusX1 = np.empty(1000)
    targetRadiusX2 = np.empty(1000)
    targetRadiusX3 = np.empty(1000)
    # condition = np.empty(1000)
    # target shape variable (power --> (x0-t0)^p)
    PX1 = np.empty(1000)
    PX2 = np.empty(1000)
    PX3 = np.empty(1000)
    # target absorbed power
    targetAbsorb = np.zeros(1000)
    ## ray variables
    # irradiance
    irrad = np.empty(1000000)
    irrado = np.empty(1000000)
    # ray velocity
    vLight = np.empty([1000000, 3])
    # ray position @ time t
    rLightT = np.empty([1000000, 3])
    # ray position @ time t+dt
    rLightTPDT = np.empty([1000000, 3])
    # ray power
    powerX1 = np.empty(1000000)
    powerX2 = np.empty(1000000)
    powerX3 = np.empty(1000000)
    # wavelength
    NTIWVariable = np.empty(1000000)
    NTIVariable = np.empty(1000000)
    wavelength = np.empty(1000000)
    # misc
    tFlag = np.zeros(1000000) # flag for if ray hits target
    iKey = np.zeros(1000000) # key to designate which target ray hits

    #################### Initializing simulation variables #####################

    # calculate number of rays per face
    nRaysPlusX1 = (rayDensity * forPlusX1SourcetubeX2 * forPlusX1SourcetubeX3).astype(int)
    nRaysMinusX1 = (rayDensity * forMinusX1SourcetubeX2 * forMinusX1SourcetubeX3).astype(int)
    nRaysPlusX2 = (rayDensity * forPlusX2SourcetubeX1 * forPlusX2SourcetubeX3).astype(int)
    nRaysMinusX2 = (rayDensity * forMinusX2SourcetubeX1 * forMinusX2SourcetubeX3).astype(int)
    nRaysPlusX3 = (rayDensity * forPlusX3SourcetubeX1 * forPlusX3SourcetubeX2).astype(int)
    nRaysMinusX3 = (rayDensity * forMinusX3SourcetubeX1 * forMinusX3SourcetubeX2).astype(int)
    # calculate total number of rays
    nRays = (nRaysPlusX1 + nRaysMinusX1 + nRaysPlusX2 + nRaysMinusX2 + \
        nRaysPlusX3 + nRaysMinusX3).astype(int)

    ####################### Targets Variables ###################################

    tarDeltaX1 = 5.0 * deltaX1 # target length (m)
    tarDeltaX2 = 5.0 * deltaX2 # target width (m)
    tarDeltaX3 = 5.0 * deltaX3 # target height (m)

    ####################### Creating Targets ####################################
    '''Create targets depending on rack orientation'''

    # print that targets are being created
    print('Creating targets...')

    count = 0 # initialize counter for target number

    # Target 1 
    if rackX2X3 == True: # if rack is in x2-x3 plane
        I1 = 0 # initialize x1 index
        X1 = tarDeltaX1 * I1 # initialize x1 position
        for I2 in range(-nTargetsX2, nTargetsX2 + 1): # loop through x2 direction
            X2 = tarDeltaX2 * I2 # calculate x2 position
            for I3 in range(-nTargetsX3, nTargetsX3 + 1): # loop through x3 direction
                X3 = tarDeltaX3 * I3 # calculate x3 position
                if (I2 == 0) or (I3 == 0): # if at center of rack
                    pass # break out of loop
                else: # if at center of rack
                    # create target
                    targetX1[count] = X1
                    targetX2[count] = X2
                    targetX3[count] = X3
                    targetRadiusX1[count] = 2.75 * deltaX1
                    targetRadiusX2[count] = 2.75 * deltaX2
                    targetRadiusX3[count] = 2.75 * deltaX3
                    PX1[count] = 6.0
                    PX2[count] = 6.0
                    PX3[count] = 6.0
                    targetAbsorb[count] = 0.0
                    # increment target number
                    count += 1 
    # Target 2
    if rackX1X3 == True: # if rack is in x1-x3 plane
        for I1 in range(-nTargetsX1, nTargetsX1 + 1): # loop through x1 direction
            X1 = tarDeltaX1 * I1 # calculate x1 position
            I2 = 0 # initialize x2 index
            X2 = tarDeltaX2 * I2 # initialize x2 position
            for I3 in range(-nTargetsX3, nTargetsX3 + 1): # loop through x3 direction
                X3 = tarDeltaX3 * I3 # calculate x3 position
                if (I1 == 0) or (I3 == 0): # if at center of rack
                    pass # break out of loop
                else: # if at center of rack
                    # create target
                    targetX1[count] = X1
                    targetX2[count] = X2
                    targetX3[count] = X3
                    targetRadiusX1[count] = 2.75 * deltaX1
                    targetRadiusX2[count] = 2.75 * deltaX2
                    targetRadiusX3[count] = 2.75 * deltaX3
                    PX1[count] = 2.0
                    PX2[count] = 2.0
                    PX3[count] = 2.0
                    targetAbsorb[count] = 0.0
                    # increment target number
                    count += 1
    # Target 3
    if rackX1X2 == True: # if rack is in x1-x2 plane
        for I1 in range(-nTargetsX1, nTargetsX1 + 1): # loop through x1 direction
            X1 = tarDeltaX1 * I1 # calculate x1 position
            for I2 in range(-nTargetsX2, nTargetsX2 + 1): # loop through x2 direction
                X2 = tarDeltaX2 * I2 # calculate x2 position
                I3 = 0 # initialize x3 index
                X3 = tarDeltaX3 * I3 # initialize x3 position
                if (I1 == 0) or (I2 == 0): # if at center of rack
                    pass # break out of loop
                else: # if at center of rack
                    # create target
                    targetX1[count] = X1
                    targetX2[count] = X2
                    targetX3[count] = X3
                    targetRadiusX1[count] = 2.75 * deltaX1
                    targetRadiusX2[count] = 2.75 * deltaX2
                    targetRadiusX3[count] = 2.75 * deltaX3
                    PX1[count] = 2.0
                    PX2[count] = 2.0
                    PX3[count] = 2.0
                    targetAbsorb[count] = 0.0
                    # increment target number
                    count += 1
    # Target 4
    if rackX1 == True: # if rack is in x1 direction
        for I1 in range(-nTargetsX1, nTargetsX1 + 1): # loop through x1 direction
            X1 = tarDeltaX1 * I1 # calculate x1 position
            I2 = 0 # initialize x2 index
            X2 = tarDeltaX2 * I2 # initialize x2 position
            I3 = 0 # initialize x3 index
            X3 = tarDeltaX3 * I3 # initialize x3 position
            # create target
            targetX1[count] = X1
            targetX2[count] = X2
            targetX3[count] = X3
            targetRadiusX1[count] = 2.75 * deltaX1
            targetRadiusX2[count] = 2.75 * deltaX2
            targetRadiusX3[count] = 2.75 * deltaX3
            PX1[count] = 2.0
            PX2[count] = 2.0
            PX3[count] = 2.0
            targetAbsorb[count] = 0.0
            # increment target number
            count += 1
    # Target 5
    if rackX2 == True: # if rack is in x2 direction
        I1 = 0 # initialize x1 index
        X1 = tarDeltaX1 * I1 # initialize x1 position
        for I2 in range(-nTargetsX2, nTargetsX2 + 1): # loop through x2 direction
            X2 = tarDeltaX2 * I2 # calculate x2 position
            I3 = 0 # initialize x3 index
            X3 = tarDeltaX3 * I3 # initialize x3 position
            # create target
            targetX1[count] = X1
            targetX2[count] = X2
            targetX3[count] = X3
            targetRadiusX1[count] = 2.75 * deltaX1
            targetRadiusX2[count] = 2.75 * deltaX2
            targetRadiusX3[count] = 2.75 * deltaX3
            PX1[count] = 2.0
            PX2[count] = 2.0
            PX3[count] = 2.0
            targetAbsorb[count] = 0.0
            # increment target number
            count += 1
    # Target 6
    if rackX3 == True: # if rack is in x3 direction
        I1 = 0 # initialize x1 index
        X1 = tarDeltaX1 * I1 # initialize x1 position
        I2 = 0 # initialize x2 index
        X2 = tarDeltaX2 * I2 # initialize x2 position
        for I3 in range(-nTargetsX3, nTargetsX3 + 1): # loop through x3 direction
            X3 = tarDeltaX3 * I3 # calculate x3 position
            # create target
            targetX1[count] = X1
            targetX2[count] = X2
            targetX3[count] = X3
            targetRadiusX1[count] = 2.75 * deltaX1
            targetRadiusX2[count] = 2.75 * deltaX2
            targetRadiusX3[count] = 2.75 * deltaX3
            PX1[count] = 2.0
            PX2[count] = 2.0
            PX3[count] = 2.0
            targetAbsorb[count] = 0.0
            # increment target number
            count += 1

    # sum up total number of targets
    nTargets = count

    # reduce array sizes to only include targets that were created
    targetX1 = targetX1[:nTargets]
    targetX2 = targetX2[:nTargets]
    targetX3 = targetX3[:nTargets]
    targetRadiusX1 = targetRadiusX1[:nTargets]
    targetRadiusX2 = targetRadiusX2[:nTargets]
    targetRadiusX3 = targetRadiusX3[:nTargets]
    PX1 = PX1[:nTargets]
    PX2 = PX2[:nTargets]
    PX3 = PX3[:nTargets]
    targetAbsorb = targetAbsorb[:nTargets]

    # print that targets have been created
    print('Targets created.')

    ############################## Defining reflection rules ##############################

    # wall locations
    wallX1MinusCutoff = -3 * (nVoxelX1 + 1) * deltaX1
    wallX1PlusCutoff = 3 * (nVoxelX1 + 1) * deltaX1
    wallX2MinusCutoff = -(nVoxelX2 + 1) * deltaX2
    wallX2PlusCutoff = (nVoxelX2 + 1) * deltaX2
    wallX3MinusCutoff = -(nVoxelX3 + 1) * deltaX3
    wallX3PlusCutoff = (nVoxelX3 + 1) * deltaX3

    # wavelength dependency
    NTIWX1Minus = 10.0
    NTIWX1Plus = 4.0
    NTIWX2Minus = 1.5
    NTIWX2Plus = 3.0
    NTIWX3Minus = 9.0
    NTIWX3Plus = 2.0
    MAGRATW = 1.0
    
    NTIX1Minus = 2.0
    NTIX1Plus = 3.0
    NTIX2Minus = 4.0
    NTIX2Plus = 5.0
    NTIX3Minus = 6.0
    NTIX3Plus = 1.3
    MAGRAT = 1.0

    wavelengthX1Minus = 10.0
    wavelengthX1Plus = 4.0
    wavelengthX2Minus = 1.5
    wavelengthX2Plus = 3.0
    wavelengthX3Minus = 9.0
    wavelengthX3Plus = 2.0

    # power dependency
    totalFacePowerX1Minus = magFacePowerX1Minus * powerMagnitude
    totalFacePowerX1Plus = magFacePowerX1Plus * powerMagnitude
    totalFacePowerX2Minus = magFacePowerX2Minus * powerMagnitude
    totalFacePowerX2Plus = magFacePowerX2Plus * powerMagnitude
    totalFacePowerX3Minus = magFacePowerX3Minus * powerMagnitude
    totalFacePowerX3Plus = magFacePowerX3Plus * powerMagnitude

    irradoX1Minus = totalFacePowerX1Minus / nRaysMinusX1
    irradoX1Plus = totalFacePowerX1Plus / nRaysPlusX1
    irradoX2Minus = totalFacePowerX2Minus / nRaysMinusX2
    irradoX2Plus = totalFacePowerX2Plus / nRaysPlusX2
    irradoX3Minus = totalFacePowerX3Minus / nRaysMinusX3
    irradoX3Plus = totalFacePowerX3Plus / nRaysPlusX3

    # misc
    ref = 0.0
    radius = 0.0001

    ############################## Define time step size ##############################

    dt = stepFactor * ((deltaX1 + deltaX2 + deltaX3) / 3) / c # time step size
    timeLimit = timeMultiplier * ((2 * nVoxelX1 * deltaX1 + 2 * nVoxelX2 * deltaX2 + 2 * nVoxelX3 * deltaX3) / 3) / c # time limit
    deltaMovieFrame = timeLimit / movieFrames # time between movie frames

    ################################# Generate rays #####################################

    # print that rays are being generated
    print('Generating rays...')

    rayCount = 0 # initialize ray counter

    # -x1 direction
    if X1MinusOn == 1: # if -x1 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[:nRaysMinusX1] = irradoX1Minus # set the first nRaysMinusX1 elements to irradoX1Minus
        irrado[:nRaysMinusX1] = irradoX1Minus # set the first nRaysMinusX1 elements to irradoX1Minus
        # randomly generate ray position depending on forMinusX1Sourcetube
        rLightT[:nRaysMinusX1,0] = wallX1PlusCutoff # initialize x1 position at the wall
        rLightT[:nRaysMinusX1,1] = randseed.uniform(-forMinusX1SourcetubeX2, forMinusX1SourcetubeX2, nRaysMinusX1) # generate random x2 position
        rLightT[:nRaysMinusX1,2] = randseed.uniform(-forMinusX1SourcetubeX3, forMinusX1SourcetubeX3, nRaysMinusX1) # generate random x3 position
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A1, A1, nRaysMinusX1) - 1.0 # generate random x1 norm
        normX2 = randseed.uniform(-A2, A2, nRaysMinusX1) # generate random x2 norm
        normX3 = randseed.uniform(-A3, A3, nRaysMinusX1) # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[:nRaysMinusX1,0] = normX1 * c # set x1 velocity
        vLight[:nRaysMinusX1,1] = normX2 * c # set x2 velocity
        vLight[:nRaysMinusX1,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[:nRaysMinusX1] = irrad[:nRaysMinusX1] * normX1 # set x1 power
        powerX2[:nRaysMinusX1] = irrad[:nRaysMinusX1] * normX2 # set x2 power
        powerX3[:nRaysMinusX1] = irrad[:nRaysMinusX1] * normX3 # set x3 power
        # wavelength
        NTIVariable[:nRaysMinusX1] = NTIX1Minus # set NTI variable
        NTIWVariable[:nRaysMinusX1] = NTIWX1Minus # set wavelength variable
        wavelength[:nRaysMinusX1] = wavelengthX1Minus # set wavelength
        # update ray count
        rayCount += nRaysMinusX1     
    # +x1 direction
    if X1PlusOn == 1: # if +x1 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[rayCount:rayCount + nRaysPlusX1] = irradoX1Plus # set the first nRaysPlusX1 elements to irradoX1Plus
        irrado[rayCount:rayCount + nRaysPlusX1] = irradoX1Plus # set the first nRaysPlusX1 elements to irradoX1Plus
        # randomly generate ray position depending on forPlusX1Sourcetube
        rLightT[rayCount:rayCount + nRaysPlusX1,0] = wallX1MinusCutoff # initialize x1 position at the wall
        rLightT[rayCount:rayCount + nRaysPlusX1,1] = randseed.uniform(-forPlusX1SourcetubeX2, forPlusX1SourcetubeX2, nRaysPlusX1) # generate random x2 position
        rLightT[rayCount:rayCount + nRaysPlusX1,2] = randseed.uniform(-forPlusX1SourcetubeX3, forPlusX1SourcetubeX3, nRaysPlusX1) # generate random x3 position
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A4, A4, nRaysPlusX1) + 1.0 # generate random x1 norm
        normX2 = randseed.uniform(-A5, A5, nRaysPlusX1) # generate random x2 norm
        normX3 = randseed.uniform(-A6, A6, nRaysPlusX1) # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[rayCount:rayCount + nRaysPlusX1,0] = normX1 * c # set x1 velocity
        vLight[rayCount:rayCount + nRaysPlusX1,1] = normX2 * c # set x2 velocity
        vLight[rayCount:rayCount + nRaysPlusX1,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[rayCount:rayCount + nRaysPlusX1] = irrad[rayCount:rayCount + nRaysPlusX1] * normX1 # set x1 power
        powerX2[rayCount:rayCount + nRaysPlusX1] = irrad[rayCount:rayCount + nRaysPlusX1] * normX2 # set x2 power
        powerX3[rayCount:rayCount + nRaysPlusX1] = irrad[rayCount:rayCount + nRaysPlusX1] * normX3 # set x3 power
        # wavelength
        NTIVariable[rayCount:rayCount + nRaysPlusX1] = NTIX1Plus # set NTI variable
        NTIWVariable[rayCount:rayCount + nRaysPlusX1] = NTIWX1Plus # set wavelength variable
        wavelength[rayCount:rayCount + nRaysPlusX1] = wavelengthX1Plus # set wavelength
        # update ray count
        rayCount += nRaysPlusX1
    # -x2 direction
    if X2MinusOn == 1: # if -x2 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[rayCount:rayCount + nRaysMinusX2] = irradoX2Minus # set the first nRaysMinusX2 elements to irradoX2Minus
        irrado[rayCount:rayCount + nRaysMinusX2] = irradoX2Minus # set the first nRaysMinusX2 elements to irradoX2Minus
        # randomly generate ray position depending on forMinusX2Sourcetube
        rLightT[rayCount:rayCount + nRaysMinusX2,0] = randseed.uniform(-forMinusX2SourcetubeX1, forMinusX2SourcetubeX1, nRaysMinusX2) # generate random x1 position
        rLightT[rayCount:rayCount + nRaysMinusX2,1] = wallX2PlusCutoff # initialize x2 position at the wall
        rLightT[rayCount:rayCount + nRaysMinusX2,2] = randseed.uniform(-forMinusX2SourcetubeX3, forMinusX2SourcetubeX3, nRaysMinusX2) # generate random x3 position
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A7, A7, nRaysMinusX2) # generate random x1 norm
        normX2 = randseed.uniform(-A8, A8, nRaysMinusX2) - 1.0 # generate random x2 norm
        normX3 = randseed.uniform(-A9, A9, nRaysMinusX2) # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[rayCount:rayCount + nRaysMinusX2,0] = normX1 * c # set x1 velocity
        vLight[rayCount:rayCount + nRaysMinusX2,1] = normX2 * c # set x2 velocity
        vLight[rayCount:rayCount + nRaysMinusX2,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[rayCount:rayCount + nRaysMinusX2] = irrad[rayCount:rayCount + nRaysMinusX2] * normX1 # set x1 power
        powerX2[rayCount:rayCount + nRaysMinusX2] = irrad[rayCount:rayCount + nRaysMinusX2] * normX2 # set x2 power
        powerX3[rayCount:rayCount + nRaysMinusX2] = irrad[rayCount:rayCount + nRaysMinusX2] * normX3 # set x3 power
        # wavelength
        NTIVariable[rayCount:rayCount + nRaysMinusX2] = NTIX2Minus # set NTI variable
        NTIWVariable[rayCount:rayCount + nRaysMinusX2] = NTIWX2Minus # set wavelength variable
        wavelength[rayCount:rayCount + nRaysMinusX2] = wavelengthX2Minus # set wavelength
        # update ray count
        rayCount += nRaysMinusX2
    # +x2 direction
    if X2PlusOn == 1: # if +x2 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[rayCount:rayCount + nRaysPlusX2] = irradoX2Plus # set the first nRaysPlusX2 elements to irradoX2Plus
        irrado[rayCount:rayCount + nRaysPlusX2] = irradoX2Plus # set the first nRaysPlusX2 elements to irradoX2Plus
        # randomly generate ray position depending on forPlusX2Sourcetube
        rLightT[rayCount:rayCount + nRaysPlusX2,0] = randseed.uniform(-forPlusX2SourcetubeX1, forPlusX2SourcetubeX1, nRaysPlusX2) # generate random x1 position
        rLightT[rayCount:rayCount + nRaysPlusX2,1] = wallX2MinusCutoff # initialize x2 position at the wall
        rLightT[rayCount:rayCount + nRaysPlusX2,2] = randseed.uniform(-forPlusX2SourcetubeX3, forPlusX2SourcetubeX3, nRaysPlusX2) # generate random x3 position
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A10, A10, nRaysPlusX2) # generate random x1 norm
        normX2 = randseed.uniform(-A11, A11, nRaysPlusX2) + 1.0 # generate random x2 norm
        normX3 = randseed.uniform(-A12, A12, nRaysPlusX2) # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[rayCount:rayCount + nRaysPlusX2,0] = normX1 * c # set x1 velocity
        vLight[rayCount:rayCount + nRaysPlusX2,1] = normX2 * c # set x2 velocity
        vLight[rayCount:rayCount + nRaysPlusX2,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[rayCount:rayCount + nRaysPlusX2] = irrad[rayCount:rayCount + nRaysPlusX2] * normX1 # set x1 power
        powerX2[rayCount:rayCount + nRaysPlusX2] = irrad[rayCount:rayCount + nRaysPlusX2] * normX2 # set x2 power
        powerX3[rayCount:rayCount + nRaysPlusX2] = irrad[rayCount:rayCount + nRaysPlusX2] * normX3 # set x3 power
        # wavelength
        NTIVariable[rayCount:rayCount + nRaysPlusX2] = NTIX2Plus # set NTI variable
        NTIWVariable[rayCount:rayCount + nRaysPlusX2] = NTIWX2Plus # set wavelength variable
        wavelength[rayCount:rayCount + nRaysPlusX2] = wavelengthX2Plus # set wavelength
        # update ray count
        rayCount += nRaysPlusX2
    # -x3 direction
    if X3MinusOn == 1: # if -x3 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[rayCount:rayCount + nRaysMinusX3] = irradoX3Minus # set the first nRaysMinusX3 elements to irradoX3Minus
        irrado[rayCount:rayCount + nRaysMinusX3] = irradoX3Minus # set the first nRaysMinusX3 elements to irradoX3Minus
        # randomly generate ray position depending on forMinusX3Sourcetube
        rLightT[rayCount:rayCount + nRaysMinusX3,0] = randseed.uniform(-forMinusX3SourcetubeX1, forMinusX3SourcetubeX1, nRaysMinusX3) # generate random x1 position
        rLightT[rayCount:rayCount + nRaysMinusX3,1] = randseed.uniform(-forMinusX3SourcetubeX2, forMinusX3SourcetubeX2, nRaysMinusX3) # generate random x2 position
        rLightT[rayCount:rayCount + nRaysMinusX3,2] = wallX3PlusCutoff # initialize x3 position at the wall
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A13, A13, nRaysMinusX3) # generate random x1 norm
        normX2 = randseed.uniform(-A14, A14, nRaysMinusX3) # generate random x2 norm
        normX3 = randseed.uniform(-A15, A15, nRaysMinusX3) - 1.0 # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[rayCount:rayCount + nRaysMinusX3,0] = normX1 * c # set x1 velocity
        vLight[rayCount:rayCount + nRaysMinusX3,1] = normX2 * c # set x2 velocity
        vLight[rayCount:rayCount + nRaysMinusX3,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[rayCount:rayCount + nRaysMinusX3] = irrad[rayCount:rayCount + nRaysMinusX3] * normX1 # set x1 power
        powerX2[rayCount:rayCount + nRaysMinusX3] = irrad[rayCount:rayCount + nRaysMinusX3] * normX2 # set x2 power
        powerX3[rayCount:rayCount + nRaysMinusX3] = irrad[rayCount:rayCount + nRaysMinusX3] * normX3 # set x3 power
        # wavelength
        NTIVariable[rayCount:rayCount + nRaysMinusX3] = NTIX3Minus # set NTI variable
        NTIWVariable[rayCount:rayCount + nRaysMinusX3] = NTIWX3Minus # set wavelength variable
        wavelength[rayCount:rayCount + nRaysMinusX3] = wavelengthX3Minus # set wavelength
        # update ray count
        rayCount += nRaysMinusX3
    # +x3 direction
    if X3PlusOn == 1: # if +x3 direction has rays
        ## initialize variables with vector operations
        # irradiance
        irrad[rayCount:rayCount + nRaysPlusX3] = irradoX3Plus # set the first nRaysPlusX3 elements to irradoX3Plus
        irrado[rayCount:rayCount + nRaysPlusX3] = irradoX3Plus # set the first nRaysPlusX3 elements to irradoX3Plus
        # randomly generate ray position depending on forPlusX3Sourcetube
        rLightT[rayCount:rayCount + nRaysPlusX3,0] = randseed.uniform(-forPlusX3SourcetubeX1, forPlusX3SourcetubeX1, nRaysPlusX3) # generate random x1 position
        rLightT[rayCount:rayCount + nRaysPlusX3,1] = randseed.uniform(-forPlusX3SourcetubeX2, forPlusX3SourcetubeX2, nRaysPlusX3) # generate random x2 position
        rLightT[rayCount:rayCount + nRaysPlusX3,2] = wallX3MinusCutoff # initialize x3 position at the wall
        # randomly generate norms to later be normalized
        normX1 = randseed.uniform(-A16, A16, nRaysPlusX3) # generate random x1 norm
        normX2 = randseed.uniform(-A17, A17, nRaysPlusX3) # generate random x2 norm
        normX3 = randseed.uniform(-A18, A18, nRaysPlusX3) + 1.0 # generate random x3 norm
        # normalize norms
        normX1 = normX1 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x1 norm
        normX2 = normX2 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x2 norm
        normX3 = normX3 / np.sqrt(normX1**2 + normX2**2 + normX3**2) # normalize x3 norm
        # ray velocity
        vLight[rayCount:rayCount + nRaysPlusX3,0] = normX1 * c # set x1 velocity
        vLight[rayCount:rayCount + nRaysPlusX3,1] = normX2 * c # set x2 velocity
        vLight[rayCount:rayCount + nRaysPlusX3,2] = normX3 * c # set x3 velocity
        # ray power
        powerX1[rayCount:rayCount + nRaysPlusX3] = irrad[rayCount:rayCount + nRaysPlusX3] * normX1 # set x1 power
        powerX2[rayCount:rayCount + nRaysPlusX3] = irrad[rayCount:rayCount + nRaysPlusX3] * normX2 # set x2 power
        powerX3[rayCount:rayCount + nRaysPlusX3] = irrad[rayCount:rayCount + nRaysPlusX3] * normX3 # set x3 power
        # wavelength
        NTIVariable[rayCount:rayCount + nRaysPlusX3] = NTIX3Plus # set NTI variable
        NTIWVariable[rayCount:rayCount + nRaysPlusX3] = NTIWX3Plus # set wavelength variable
        wavelength[rayCount:rayCount + nRaysPlusX3] = wavelengthX3Plus # set wavelength
        # update ray count
        rayCount += nRaysPlusX3

    # reduce array sizes to only include rays that have been created
    irrad = irrad[:rayCount] # reduce irrad array size
    irrado = irrado[:rayCount] # reduce irrado array size
    rLightT = rLightT[:rayCount,:] # reduce rLightT array size
    rLightTPDT = rLightTPDT[:rayCount,:] # reduce rLightTPDT array size
    vLight = vLight[:rayCount,:] # reduce vLight array size
    powerX1 = powerX1[:rayCount] # reduce powerX1 array size
    powerX2 = powerX2[:rayCount] # reduce powerX2 array size
    powerX3 = powerX3[:rayCount] # reduce powerX3 array size
    NTIVariable = NTIVariable[:rayCount] # reduce NTIVariable array size
    NTIWVariable = NTIWVariable[:rayCount] # reduce NTIWVariable array size
    wavelength = wavelength[:rayCount] # reduce wavelength array size
    
    # print that the rays have been created
    print('Rays created.')

    ############################ Propagate rays #########################################

    # print line to separate data
    print('----------------------------------------')
    print('Propagating rays...')

    # initialize time loop
    time = 0.0 # initialize time
    while time < timeLimit: # while time is less than timeLimit, propagate rays

        # print time every 10 timesteps
        # if time % 10 == 0:
        print('Time:',  time, 's.')

        ############################## no reflections ########################################

        if reflections == False: # with no ray reflections
            ## check for ray - target collisions by iterating through each target
            for i in range(nTargets): # iterate through each target
                # check if rays are within the target
                condition = (abs(rLightT[:,0] - targetX1[i]) / targetRadiusX1[i])**PX1[i] + \
                    (abs(rLightT[:,1] - targetX2[i]) / targetRadiusX2[i])**PX2[i] + \
                    (abs(rLightT[:,2] - targetX3[i]) / targetRadiusX3[i])**PX3[i] <= 1.0
                # flag rays that are within the target
                tFlag[np.where(condition)] = 1.0
                # update position of rays
                rLightTPDT = rLightT + vLight * dt

        else: # with ray reflections

            ###################### reset variables  ##############################

            # reset target flags
            # tFlag = np.zeros(rayCount) # initialize target flag array
            iKey = np.zeros(rayCount, dtype=int) # initialize target key array

            # reset norms
            NX1 = np.zeros(rayCount) # initialize x1 norm array
            NX2 = np.zeros(rayCount) # initialize x2 norm array
            NX3 = np.zeros(rayCount) # initialize x3 norm array

            ###################### check for ray - target collisions ##############################

            # check for ray - target collisions by iterating through each target
            for i in range(nTargets): # iterate through each target
                # check if rays are within the target
                condition = (abs(rLightT[:,0] - targetX1[i]) / targetRadiusX1[i])**PX1[i] + \
                    (abs(rLightT[:,1] - targetX2[i]) / targetRadiusX2[i])**PX2[i] + \
                    (abs(rLightT[:,2] - targetX3[i]) / targetRadiusX3[i])**PX3[i] <= 1.0
                # flag rays that are within the target
                tFlag[np.where(condition)[0]] = 1.0
                # flag which target the ray is in
                iKey[np.where(condition)[0]] = i + 1

            # get indices of rays that are within a target
            refTargetInd = np.where(iKey != 0)[0]

            ############################## check for ray - wall collisions ########################

            # flag rays that are within a wall (for each statement, make sure that the ray is not already flagged)
            wX1PlusFlag = rLightT[:,0] >= wallX1PlusCutoff # -x1 wall
            wX1MinusFlag = rLightT[:,0] <= wallX1MinusCutoff # +x1 wall
            wX2PlusFlag = (rLightT[:,1] >= wallX2PlusCutoff) * (wX1MinusFlag == False) * (wX1PlusFlag == False) # -x2 wall
            wX2MinusFlag = (rLightT[:,1] <= wallX2MinusCutoff) * (wX1MinusFlag == False) * (wX1PlusFlag == False) # +x2 wall
            wX3PlusFlag = (rLightT[:,2] >= wallX3PlusCutoff) * (wX1MinusFlag == False) * (wX1PlusFlag == False) * \
                (wX2MinusFlag == False) * (wX2PlusFlag == False) # -x3 wall
            wX3MinusFlag = (rLightT[:,2] <= wallX3MinusCutoff) * (wX1MinusFlag == False) * (wX1PlusFlag == False) * \
                (wX2MinusFlag == False) * (wX2PlusFlag == False) # +x3 wall
            
            # a general flag for rays that are within a wall
            wFlag = wX1MinusFlag + wX1PlusFlag + wX2MinusFlag + wX2PlusFlag + wX3MinusFlag + wX3PlusFlag

            # calculate updated norms for rays that are within a wall
            NX1 = wX1PlusFlag * -1.0 + wX1MinusFlag * 1.0 + wX2PlusFlag * 0.0 + wX2MinusFlag * 0.0 + \
                wX3PlusFlag * 0.0 + wX3MinusFlag * 0.0
            NX2 = wX1PlusFlag * 0.0 + wX1MinusFlag * 0.0 + wX2PlusFlag * -1.0 + wX2MinusFlag * 1.0 + \
                wX3PlusFlag * 0.0 + wX3MinusFlag * 0.0
            NX3 = wX1PlusFlag * 0.0 + wX1MinusFlag * 0.0 + wX2PlusFlag * 0.0 + wX2MinusFlag * 0.0 + \
                wX3PlusFlag * -1.0 + wX3MinusFlag * 1.0
            
            ############################## if rays are within a target ###########################

            if refTargetInd.size > 0: # if there are rays within a target

                # get number of rays within a target
                nRefTarget = refTargetInd.size

                gradF = np.zeros((nRefTarget,3)) # initialize gradient of F

                # target indices per ray that is within a target
                tpr = iKey[refTargetInd] - 1 # target per ray (tpr)

                # calculate distance to target
                xDist = rLightT[refTargetInd,0] - targetX1[tpr] # x distance
                yDist = rLightT[refTargetInd,1] - targetX2[tpr] # y distance
                zDist = rLightT[refTargetInd,2] - targetX3[tpr] # z distance

                # target indices for rays that satisfy abs(xDist) > 0 / abs(yDist) > 0 / abs(zDist) > 0

                xDistInd = np.where(abs(xDist) > 0.00000001)[0]
                yDistInd = np.where(abs(yDist) > 0.00000001)[0]
                zDistInd = np.where(abs(zDist) > 0.00000001)[0]

                indGradF1 = tpr[xDistInd]
                indGradF2 = tpr[yDistInd]
                indGradF3 = tpr[zDistInd]
                
                # computing normals
                if indGradF1.size != 0: # if there are rays that satisfy abs(xDist) > 0
                    term1 = PX1[indGradF1] * ((abs(xDist[xDistInd]) / targetRadiusX1[indGradF1])**(PX1[indGradF1] - 1.0))
                    term2 = xDist[xDistInd] / (targetRadiusX1[indGradF1] * abs(xDist[xDistInd]))
                    gradF[xDistInd,0] = term1 * term2 # x1 component of gradient of F
                
                if indGradF2.size != 0: # if there are rays that satisfy abs(yDist) > 0
                    term1 = PX2[indGradF2] * ((abs(yDist[yDistInd]) / targetRadiusX2[indGradF2])**(PX2[indGradF2] - 1.0))
                    term2 = yDist[yDistInd] / (targetRadiusX2[indGradF2] * abs(yDist[yDistInd]))
                    gradF[yDistInd,1] = term1 * term2 # x2 component of gradient of F

                if indGradF3.size != 0: # if there are rays that satisfy abs(zDist) > 0
                    term1 = PX3[indGradF3] * ((abs(zDist[zDistInd]) / targetRadiusX3[indGradF3])**(PX3[indGradF3] - 1.0))
                    term2 = zDist[zDistInd] / (targetRadiusX3[indGradF3] * abs(zDist[zDistInd]))
                    gradF[zDistInd,2] = term1 * term2 # x3 component of gradient of F

                # calculate the normal vector
                gNorm = np.sqrt(gradF[:,0]**2 + gradF[:,1]**2 + gradF[:,2]**2) # gradient norm

                gradF[:,0] = gradF[:,0] / gNorm # x1 component of normal vector
                gradF[:,1] = gradF[:,1] / gNorm # x2 component of normal vector
                gradF[:,2] = gradF[:,2] / gNorm # x3 component of normal vector

                NX1[refTargetInd] = gradF[:,0] # x1 component of normal vector
                NX2[refTargetInd] = gradF[:,1] # x2 component of normal vector
                NX3[refTargetInd] = gradF[:,2] # x3 component of normal vector
            
            # determining the normal projection and subtracting it out
            vNorm = np.sqrt(vLight[:,0]**2 + vLight[:,1]**2 + vLight[:,2]**2) # velocity norm
            vN = vLight[:,0] * NX1 + vLight[:,1] * NX2 + vLight[:,2] * NX3 # normal velocity

            ############################## TARGET REFLECTIONS ##############################
            ################################################################################

            ## if vN is or less than zero, then the ray is being reflected
            refIndT = (vN <= 0.0) * (iKey != 0.0) # get the indices of the rays that are being reflected and are moving towards a target
            # get the indices of the rays that are being reflected and are moving towards a target
            refIndT = np.where(refIndT)[0]

            if refIndT.size > 0: # if there are rays that are moving towards a target

                # calculate the incident angle
                thetaI = np.arccos(abs(vN[refIndT]) / vNorm[refIndT]) # incident angle

                # reflect rays
                vLight[refIndT,0] = vLight[refIndT,0] - 2.0 * NX1[refIndT] * vN[refIndT] # x1 velocity
                vLight[refIndT,1] = vLight[refIndT,1] - 2.0 * NX2[refIndT] * vN[refIndT] # x2 velocity
                vLight[refIndT,2] = vLight[refIndT,2] - 2.0 * NX3[refIndT] * vN[refIndT] # x3 velocity

                # calculate reflectivity (IR)
                rPer = (abs(np.sin(thetaI) / NTIVariable[refIndT]) > 1.0) * 1.0 + \
                    (abs(np.sin(thetaI) / NTIVariable[refIndT]) <= 1.0) * \
                    (np.cos(thetaI) - (1.0 / MAGRAT) * (NTIVariable[refIndT]**2 - np.sin(thetaI)**2)**0.5) / \
                    (np.cos(thetaI) + (1.0 / MAGRAT) * (NTIVariable[refIndT]**2 - np.sin(thetaI)**2)**0.5)
                
                rPar = (abs(np.sin(thetaI) / NTIVariable[refIndT]) > 1.0) * 1.0 + \
                    (abs(np.sin(thetaI) / NTIVariable[refIndT]) <= 1.0) * \
                    ((1.0 / MAGRAT) * (NTIVariable[refIndT]**2) * np.cos(thetaI) - \
                    (NTIVariable[refIndT]**2 - np.sin(thetaI)**2)**0.5) / \
                    ((1.0 / MAGRAT) * (NTIVariable[refIndT]**2) * np.cos(thetaI) + \
                    (NTIVariable[refIndT]**2 - np.sin(thetaI)**2)**0.5)
                
                ############################## reflectivities ########################################

                bigRPar = rPar**2 # parallel reflectivity
                bigRPer = rPer**2 # perpendicular reflectivity
                bigR = 0.5 * (bigRPar + bigRPer) # average reflectivity

                ######################## absorption and irradiance ###################################
                
                # calculate absorption by looping through the targets
                for i in range(nTargets):
                    # sum up contributions of each ray to the target
                    targetAbsorb[i] += np.sum(irrad[refIndT]* (1 - bigR) * (iKey[refIndT] == i+1))

                irrad[refIndT] = irrad[refIndT] * bigR # irradiance

                # update ray power
                powerX1[refIndT] = vLight[refIndT,0] * irrad[refIndT] / vNorm[refIndT] # set x1 power
                powerX2[refIndT] = vLight[refIndT,1] * irrad[refIndT] / vNorm[refIndT] # set x2 power
                powerX3[refIndT] = vLight[refIndT,2] * irrad[refIndT] / vNorm[refIndT] # set x3 power

                
            ############################## WALL REFLECTIONS ########################################  
            ########################################################################################

            ## if vN is or less than zero, then the ray is being reflected
            refInd = (vN <= 0.0) * (wFlag == 1.0) # get the indices of the rays that are being reflected and are moving towards the wall
            # get the indices of the rays that are being reflected and are moving towards the wall
            refInd = np.where(refInd)[0]

            if refInd.size > 0: # if there are rays that are moving towards the wall

                # calculate the incident angle
                thetaI = np.arccos(abs(vN[refInd]) / vNorm[refInd]) # incident angle

                # reflect rays
                vLight[refInd,0] = vLight[refInd,0] - 2.0 * NX1[refInd] * vN[refInd] # x1 velocity
                vLight[refInd,1] = vLight[refInd,1] - 2.0 * NX2[refInd] * vN[refInd] # x2 velocity
                vLight[refInd,2] = vLight[refInd,2] - 2.0 * NX3[refInd] * vN[refInd] # x3 velocity

                # calculate reflectivity (IR)
                rPer = (abs(np.sin(thetaI) / NTIWVariable[refInd]) > 1.0) * 1.0 + \
                    (abs(np.sin(thetaI) / NTIWVariable[refInd]) <= 1.0) * \
                    (np.cos(thetaI) - (1.0 / MAGRATW) * (NTIWVariable[refInd]**2 - np.sin(thetaI)**2)**0.5) / \
                    (np.cos(thetaI) + (1.0 / MAGRATW) * (NTIWVariable[refInd]**2 - np.sin(thetaI)**2)**0.5)
                
                rPar = (abs(np.sin(thetaI) / NTIWVariable[refInd]) > 1.0) * 1.0 + \
                    (abs(np.sin(thetaI) / NTIWVariable[refInd]) <= 1.0) * \
                    ((1.0 / MAGRATW) * (NTIWVariable[refInd]**2) * np.cos(thetaI) - \
                    (NTIWVariable[refInd]**2 - np.sin(thetaI)**2)**0.5) / \
                    ((1.0 / MAGRATW) * (NTIWVariable[refInd]**2) * np.cos(thetaI) + \
                    (NTIWVariable[refInd]**2 - np.sin(thetaI)**2)**0.5)
                

                # get indices of rays that are moving towards each wall
                wX1PlusInd = np.where((vN <= 0.0) * (wX1PlusFlag== 1.0))[0] # +x1 wall
                wX1MinusInd = np.where((vN <= 0.0) * (wX1MinusFlag== 1.0))[0] # -x1 wall
                wX2PlusInd = np.where((vN <= 0.0) * (wX2PlusFlag== 1.0))[0] # +x2 wall
                wX2MinusInd = np.where((vN <= 0.0) * (wX2MinusFlag== 1.0))[0] # -x2 wall
                wX3PlusInd = np.where((vN <= 0.0) * (wX3PlusFlag== 1.0))[0] # +x3 wall
                wX3MinusInd = np.where((vN <= 0.0) * (wX3MinusFlag== 1.0))[0] # -x3 wall
                                                                   
                # reset location to start at the wall with an offset
                rLightT[wX1PlusInd,0] = wallX1PlusCutoff - 1.5 * deltaX1 
                rLightT[wX1MinusInd,0] = wallX1MinusCutoff + 1.5 * deltaX1
                rLightT[wX2PlusInd,1] = wallX2PlusCutoff - 1.5 * deltaX2
                rLightT[wX2MinusInd,1] = wallX2MinusCutoff + 1.5 * deltaX2
                rLightT[wX3PlusInd,2] = wallX3PlusCutoff - 1.5 * deltaX3
                rLightT[wX3MinusInd,2] = wallX3MinusCutoff + 1.5 * deltaX3

                ############################## reflectivities ########################################

                bigRPar = rPar**2 # parallel reflectivity
                bigRPer = rPer**2 # perpendicular reflectivity
                bigR = 0.5 * (bigRPar + bigRPer) # average reflectivity

                ######################## absorption and irradiance ###################################

                irrad[refInd] = irrad[refInd] * bigR # irradiance

                # update ray power
                powerX1[refInd] = vLight[refInd,0] * irrad[refInd] / vNorm[refInd] # set x1 power
                powerX2[refInd] = vLight[refInd,1] * irrad[refInd] / vNorm[refInd] # set x2 power
                powerX3[refInd] = vLight[refInd,2] * irrad[refInd] / vNorm[refInd] # set x3 power


            ############################## update position ########################################
            # (1. if the power is below rayTol --> trap the light since it is absorbed,
            # 2. otherwise --> update position using the velocity)
            rLightTPDT[:,0] = (irrad / irrado <= rayTol) * rLightT[:,0] + \
                (irrad / irrado > rayTol) * (rLightT[:,0] + vLight[:,0] * dt)
            rLightTPDT[:,1] = (irrad / irrado <= rayTol) * rLightT[:,1] + \
                (irrad / irrado > rayTol) * (rLightT[:,1] + vLight[:,1] * dt)
            rLightTPDT[:,2] = (irrad / irrado <= rayTol) * rLightT[:,2] + \
                (irrad / irrado > rayTol) * (rLightT[:,2] + vLight[:,2] * dt)

        ############################ Write out data file if toggled ############################

        # save data file if time == movieTime
        if time >= movieTime:

            if raySave == True: # if raySave is toggled on, save ray data file

            # The ray variables are written out in the following order: 
            # x-pos, y-pos, z-pos, radius, irradiance ratio (IR/IRo), x-vel, y-vel, z-vel, 
            # powerX1/IRo, powerX2,IRo, powerX3/IRo, 0.0, wavelength
                                
                #initialize data file
                with open('ray.dat', 'a') as p_file:
                    writer = csv.writer(p_file, delimiter='\t',quotechar = "'")  # make write variable
                    if time == 0.0: #only write title and heading at first time step
                        writer.writerow(['TITLE = ', 'Lightbox Simulation'])
                        writer.writerow(['VARIABLES = ', 'X-coord', 'Y-coord', 'Z-coord','TargetR', 'IR/IRo', 'X-vel', 'Y-vel', 'Z-vel', 'PowerX1/IRo', 'PowerX2/IRo', 'PowerX3/IRo', 'Var', 'Wavelength'])  
                    
                    #write position and radius data
                    writer.writerow(['ZONE T = \"RAY\"']) # write zone title
                    writer.writerow(['STRANDID=',2]) # write strand ID
                    writer.writerow(['SOLUTIONTIME=',time]) # write solution time
                    for i in range(0, rayCount):
                        writer.writerow([rLightTPDT[i,0], rLightTPDT[i,1], rLightTPDT[i,2], 2.0 * radius, irrad[i]/irrado[i], vLight[i,0],\
                                        vLight[i,1], vLight[i,2], powerX1[i]/irrado[i], powerX2[i]/irrado[i], powerX3[i]/irrado[i],\
                                            0.0, wavelength[i]])
                        
            # write out target data file
            if targetSave == True:
                                
                #initialize data file
                with open('target.dat', 'a') as p_file:
                    writer = csv.writer(p_file, delimiter='\t',quotechar = "'")  # make write variable
                    # if i == 0: #only write title and heading at first time step
                    writer.writerow(['TITLE = ', 'Lightbox Simulation'])
                    writer.writerow(['VARIABLES = ', 'X-coord', 'Y-coord', 'Z-coord','TargetR','Absorb-Ratio'])  

                    #write position and radius data
                    writer.writerow(['ZONE T = \"TARGET\"']) # write zone title
                    writer.writerow(['STRANDID=',1]) # write strand ID
                    writer.writerow(['SOLUTIONTIME=',time]) # write solution time
                    for i in range(0, nTargets):
                        writer.writerow([targetX1[i], targetX2[i], targetX3[i], 2.0 * targetRadiusX1[i], \
                                         targetAbsorb[i]/(totalFacePowerX1Minus + totalFacePowerX1Plus + \
                                                            totalFacePowerX2Minus + totalFacePowerX2Plus + \
                                                                totalFacePowerX3Minus + totalFacePowerX3Plus)]) # write data

                # update movieTime
                movieTime += deltaMovieFrame # update movieTime


        # update variables for next time step
        rLightT = rLightTPDT # update position

        # iterate time
        time += dt # iterate time
    
    ############################ Print important information ############################
    # print line to separate data
    print('----------------------------------------')

    # calculate number of rays hitting the target by summing the number of rays that are flagged
    rCount = np.sum(tFlag)
    rLights = rayCount # number of rays
    sum = np.sum(targetAbsorb) # sum of absorbed power
    absorbDiff = np.abs(max(targetAbsorb) - np.mean(targetAbsorb)) / max(targetAbsorb) # difference between max and min absorbed power (normalized)

    print('Rays hitting target: %d' % rCount) # print number of rays hitting target
    print('Total rays: ', rLights) # print out number of rays
    print('Fraction of rays hitting target: %f' % (rCount / rLights)) # print fraction of rays hitting target
    print('Number of Targets: %d' % nTargets) # print number of targets
    print('Fraction of power absorbed by pods: %f' % (sum / (totalFacePowerX1Minus * X1MinusOn + totalFacePowerX1Plus * X1PlusOn + \
                                                            totalFacePowerX2Minus * X2MinusOn + totalFacePowerX2Plus * X2PlusOn + \
                                                                totalFacePowerX3Minus * X3MinusOn + totalFacePowerX3Plus * X3PlusOn))) # print fraction of power absorbed by pods
    print('Power usage effectiveness by pods: %f' % ((totalFacePowerX1Minus * X1MinusOn + totalFacePowerX1Plus * X1PlusOn + \
                                                            totalFacePowerX2Minus * X2MinusOn + totalFacePowerX2Plus * X2PlusOn + \
                                                                totalFacePowerX3Minus * X3MinusOn + totalFacePowerX3Plus * X3PlusOn) / sum)) # print power usage effectiveness by pods
    print('Target absorption discrepancy: %f' % (absorbDiff)) # print target absorption discrepancy

    # CALCULATE TOTAL COST (OBJECTIVE FUNCTION)

    # calculate cost parameter associated with fraction of power absorbed by pods
    alpha = 1.0 - (sum / (totalFacePowerX1Minus * X1MinusOn + totalFacePowerX1Plus * X1PlusOn + \
                                                            totalFacePowerX2Minus * X2MinusOn + totalFacePowerX2Plus * X2PlusOn + \
                                                                totalFacePowerX3Minus * X3MinusOn + totalFacePowerX3Plus * X3PlusOn))
    
    # calculate cost parameter associated with target absorption discrepancy
    gamma = absorbDiff

    # calculate total cost
    Pi = alpha

    print('Total Cost: %f' % Pi) # print total cost

    # save outputs to dictionary
    outputs = {
        'Pi': Pi, # Cost for given input design string
        # 'alpha': alpha, # Cost parameter associated with fraction of power absorbed by pods
        # 'gamma': gamma, # Cost parameter associated with target absorption discrepancy
    }
    
    return (outputs)

# GA Function
def myGA(S,G,P,SB,Func,consts):

    # Movie settings
    movieSettings = {
        
        # toggle what to save/plot
        'targetSave': False, # save target data
        'raySave': False, # save ray data

        'movieFrames': 100, # number of frames in movie
        'movieTime': 0.0 # start time of movie [s]

}
    
    # Get number of design variables from size of search bound array
    dv = SB.shape[0] 
    
    # set number of kids (K) equal to number of parents
    K = P
    
    # Initialize all variables to be saved
    Min = np.zeros(G) # Minimum cost for each generation
    PAve = np.zeros(G) # Parent average for each generation
    Ave = np.zeros(G) # Total population average for each generation
    lamHist = np.zeros((G,dv))  # history of best design string per generation  
    
    
    Pi = np.zeros(S) # All costs in an individual generation
    
    # Initialize Lam array which will contain all of the design strings for all design variables
    Lam = np.zeros((S,dv))
    
    # Generate initial random population by looping through the number of design variables and building out the Lam array
    for i in range(dv):
        Lam[:,i] = uniform(SB[i,0],SB[i,1]-SB[i,0],size = S)
        
    # In first generation, calculate cost for all strings.
    # After, only calculate new strings since fitness for top P parents are already calculated
    # Initialize an index parameter to 0 to track which design string to start evaluating costs from
    start = 0 
    
    for i in range(G): # Loop through generations
        
        # Calculate fitness of unknown design string costs
        for j in range(start,S): # Evaluate fitness of strings
            
            # Plug in design string control variables and array of function constants
            output = Func(Lam[j,:], consts, movieSettings) # Outputs dict of function outputs
            
            Pi[j] = output['Pi'] # Extract cost from dict of outputs and assign it to cost array
            
        
        # Sort cost and design strings based on performance
        ind = np.argsort(Pi)
        Pi = np.sort(Pi)
        Lam = Lam[ind,:]
        
        # Save best design string for current generation
        lamHist[i,:] = Lam[0,:]
        
        # Generate offspring radnom parameters and indices for vectorized offspring calculation
        phi = np.random.rand(K,SB.shape[0]) # Generate random weights for offspring
        ind1 = range(0,K,2) # First set of children based on even numbered parameters
        ind2 = range(1,K,2) # Second set of children based on odd numbered parameters
        
        Parents = Lam[0:P,:] # Top P performing parents
        Children1 = phi[ind1,:]*Lam[ind1,:] + (1-phi[ind1,:])*Lam[ind2,:] # First set of children
        Children2 = phi[ind2,:]*Lam[ind2,:] + (1-phi[ind2,:])*Lam[ind1,:] # Second set of children
        
        # Initialize newPopulation array which will have S-P-K new random strings for all design variables
        newPop = np.zeros((S-P-K,dv))
        
        # Generate S - P - K new random strings by looping through the number of design variables and building out the new Population array
        for j in range(dv):
            newPop[:,j] = uniform(SB[j,0],SB[j,1]-SB[j,0],size = S-P-K)
        
         # Vertically stack parents, children, and new strings to use in next generation  
        Lam = np.vstack((Parents, Children1, Children2, newPop))  
        
        # Save minimum, parent average, and population average cost values for plotting
        Min[i] = Pi[0]
        PAve[i] = np.mean(Pi[0:P])
        Ave[i] = np.mean(Pi)
        
        # Update index parameter to P such that only new string cost values are calculated
        start = P
        
        # print line 
        print("--------------------------------------------------")
        # Print mininum value of cost for debugging (should monotonically decrease over generations)
        print("Best cost for generation " + str(i+1) + ": " + str(Min[i]))
        
    bestLam = Lam[0,:] # Extract best design string parameters afer all generations are run
    
    return(lamHist, Lam, bestLam, Pi, Min, PAve, Ave)

#%% Optimization Parameters (Test Values)

# beam spread parameters (18 total)
A1 = 0.5
A2 = 0.5
A3 = 0.5
A4 = 0.5
A5 = 0.5
A6 = 0.5
A7 = 0.5
A8 = 0.5
A9 = 0.5
A10 = 0.5
A11 = 0.5
A12 = 0.5
A13 = 0.5
A14 = 0.5
A15 = 0.5
A16 = 0.5
A17 = 0.5
A18 = 0.5

# sourcetube parameters (12 total)
forPlusX3SourcetubeX1 = 3.0
forPlusX3SourcetubeX2 = 0.5
forMinusX3SourcetubeX1 = 3.0
forMinusX3SourcetubeX2 = 0.5

forPlusX2SourcetubeX1 = 3.0
forPlusX2SourcetubeX3 = 0.5
forMinusX2SourcetubeX1 = 3.0
forMinusX2SourcetubeX3 = 0.5

forPlusX1SourcetubeX2 = 0.5
forPlusX1SourcetubeX3 = 0.5
forMinusX1SourcetubeX2 = 0.5
forMinusX1SourcetubeX3 = 0.5

# power dependency parameters (6 total)
magFacePowerX1Minus = 1e7
magFacePowerX1Plus = 1e7
magFacePowerX2Minus = 1e7
magFacePowerX2Plus = 1e7
magFacePowerX3Minus = 1e7
magFacePowerX3Plus = 1e7

# total number of design variables
numLam = 36

# Initialize design variables
Lam = np.zeros(numLam)

# Assigning design variables to associated parameters
Lam[0] = A1
Lam[1] = A2
Lam[2] = A3
Lam[3] = A4
Lam[4] = A5
Lam[5] = A6
Lam[6] = A7
Lam[7] = A8
Lam[8] = A9
Lam[9] = A10
Lam[10] = A11
Lam[11] = A12
Lam[12] = A13
Lam[13] = A14
Lam[14] = A15
Lam[15] = A16
Lam[16] = A17
Lam[17] = A18
Lam[18] = forPlusX3SourcetubeX1
Lam[19] = forPlusX3SourcetubeX2
Lam[20] = forMinusX3SourcetubeX1
Lam[21] = forMinusX3SourcetubeX2
Lam[22] = forPlusX2SourcetubeX1
Lam[23] = forPlusX2SourcetubeX3
Lam[24] = forMinusX2SourcetubeX1
Lam[25] = forMinusX2SourcetubeX3
Lam[26] = forPlusX1SourcetubeX2
Lam[27] = forPlusX1SourcetubeX3
Lam[28] = forMinusX1SourcetubeX2
Lam[29] = forMinusX1SourcetubeX3
Lam[30] = magFacePowerX1Minus
Lam[31] = magFacePowerX1Plus
Lam[32] = magFacePowerX2Minus
Lam[33] = magFacePowerX2Plus
Lam[34] = magFacePowerX3Minus
Lam[35] = magFacePowerX3Plus

# run lightbox simulation
Pi = lightboxSim(Lam, lightboxConstants, movieSettings)
# %% ################# set lower and upper bounds for design variables #################

# numLam = 36 # total number of design variables

# ## lower limits

# # beam spread parameters (18 total)
# A1Min = 0.0
# A2Min = 0.0
# A3Min = 0.0
# A4Min = 0.0
# A5Min = 0.0
# A6Min = 0.0
# A7Min = 0.0
# A8Min = 0.0
# A9Min = 0.0
# A10Min = 0.0
# A11Min = 0.0
# A12Min = 0.0
# A13Min = 0.0
# A14Min = 0.0
# A15Min = 0.0
# A16Min = 0.0
# A17Min = 0.0
# A18Min = 0.0

# # sourcetube parameters (12 total)
# forPlusX3SourcetubeX1Min = 0.0
# forPlusX3SourcetubeX2Min = 0.0
# forMinusX3SourcetubeX1Min = 0.0
# forMinusX3SourcetubeX2Min = 0.0

# forPlusX2SourcetubeX1Min = 0.0
# forPlusX2SourcetubeX3Min = 0.0
# forMinusX2SourcetubeX1Min = 0.0
# forMinusX2SourcetubeX3Min = 0.0

# forPlusX1SourcetubeX2Min = 0.0
# forPlusX1SourcetubeX3Min = 0.0
# forMinusX1SourcetubeX2Min = 0.0
# forMinusX1SourcetubeX3Min = 0.0

# # power dependency parameters (6 total)
# magFacePowerX1MinusMin = 1e6
# magFacePowerX1PlusMin = 1e6
# magFacePowerX2MinusMin = 1e6
# magFacePowerX2PlusMin = 1e6
# magFacePowerX3MinusMin = 1e6
# magFacePowerX3PlusMin = 1e6

# # stack lower limits
# lowBound = np.zeros([numLam,1])

# lowBound[0] = A1Min
# lowBound[1] = A2Min
# lowBound[2] = A3Min
# lowBound[3] = A4Min
# lowBound[4] = A5Min
# lowBound[5] = A6Min
# lowBound[6] = A7Min
# lowBound[7] = A8Min
# lowBound[8] = A9Min
# lowBound[9] = A10Min
# lowBound[10] = A11Min
# lowBound[11] = A12Min
# lowBound[12] = A13Min
# lowBound[13] = A14Min
# lowBound[14] = A15Min
# lowBound[15] = A16Min
# lowBound[16] = A17Min
# lowBound[17] = A18Min
# lowBound[18] = forPlusX3SourcetubeX1Min
# lowBound[19] = forPlusX3SourcetubeX2Min
# lowBound[20] = forMinusX3SourcetubeX1Min
# lowBound[21] = forMinusX3SourcetubeX2Min
# lowBound[22] = forPlusX2SourcetubeX1Min
# lowBound[23] = forPlusX2SourcetubeX3Min
# lowBound[24] = forMinusX2SourcetubeX1Min
# lowBound[25] = forMinusX2SourcetubeX3Min
# lowBound[26] = forPlusX1SourcetubeX2Min
# lowBound[27] = forPlusX1SourcetubeX3Min
# lowBound[28] = forMinusX1SourcetubeX2Min
# lowBound[29] = forMinusX1SourcetubeX3Min
# lowBound[30] = magFacePowerX1MinusMin
# lowBound[31] = magFacePowerX1PlusMin
# lowBound[32] = magFacePowerX2MinusMin
# lowBound[33] = magFacePowerX2PlusMin
# lowBound[34] = magFacePowerX3MinusMin
# lowBound[35] = magFacePowerX3PlusMin

# ## upper limits

# # beam spread parameters (18 total)
# A1Max = 1.0
# A2Max = 1.0
# A3Max = 1.0
# A4Max = 1.0
# A5Max = 1.0
# A6Max = 1.0
# A7Max = 1.0
# A8Max = 1.0
# A9Max = 1.0
# A10Max = 1.0
# A11Max = 1.0
# A12Max = 1.0
# A13Max = 1.0
# A14Max = 1.0
# A15Max = 1.0
# A16Max = 1.0
# A17Max = 1.0
# A18Max = 1.0

# # sourcetube parameters (12 total)
# forPlusX3SourcetubeX1Max = 3.0
# forPlusX3SourcetubeX2Max = 0.5
# forMinusX3SourcetubeX1Max = 3.0
# forMinusX3SourcetubeX2Max = 0.5

# forPlusX2SourcetubeX1Max = 3.0
# forPlusX2SourcetubeX3Max = 0.5
# forMinusX2SourcetubeX1Max = 3.0
# forMinusX2SourcetubeX3Max = 0.5

# forPlusX1SourcetubeX2Max = 0.5
# forPlusX1SourcetubeX3Max = 0.5
# forMinusX1SourcetubeX2Max = 0.5
# forMinusX1SourcetubeX3Max = 0.5

# # power dependency parameters (6 total)
# magFacePowerX1MinusMax = 1e7
# magFacePowerX1PlusMax = 1e7
# magFacePowerX2MinusMax = 1e7
# magFacePowerX2PlusMax = 1e7
# magFacePowerX3MinusMax = 1e7
# magFacePowerX3PlusMax = 1e7

# # stack upper bounds
# upBound = np.zeros([numLam,1])

# upBound[0] = A1Max
# upBound[1] = A2Max
# upBound[2] = A3Max
# upBound[3] = A4Max
# upBound[4] = A5Max
# upBound[5] = A6Max
# upBound[6] = A7Max
# upBound[7] = A8Max
# upBound[8] = A9Max
# upBound[9] = A10Max
# upBound[10] = A11Max
# upBound[11] = A12Max
# upBound[12] = A13Max
# upBound[13] = A14Max
# upBound[14] = A15Max
# upBound[15] = A16Max
# upBound[16] = A17Max
# upBound[17] = A18Max
# upBound[18] = forPlusX3SourcetubeX1Max
# upBound[19] = forPlusX3SourcetubeX2Max
# upBound[20] = forMinusX3SourcetubeX1Max
# upBound[21] = forMinusX3SourcetubeX2Max
# upBound[22] = forPlusX2SourcetubeX1Max
# upBound[23] = forPlusX2SourcetubeX3Max
# upBound[24] = forMinusX2SourcetubeX1Max
# upBound[25] = forMinusX2SourcetubeX3Max
# upBound[26] = forPlusX1SourcetubeX2Max
# upBound[27] = forPlusX1SourcetubeX3Max
# upBound[28] = forMinusX1SourcetubeX2Max
# upBound[29] = forMinusX1SourcetubeX3Max
# upBound[30] = magFacePowerX1MinusMax
# upBound[31] = magFacePowerX1PlusMax
# upBound[32] = magFacePowerX2MinusMax
# upBound[33] = magFacePowerX2PlusMax
# upBound[34] = magFacePowerX3MinusMax
# upBound[35] = magFacePowerX3PlusMax

# #%% Genetic Algorithm Parameters
# K = 4 # Strings generated by breeding
# P = 4 # Surviving strings for breeding
# S = 10 # Design strings per generation
# G = 50 # Total Generations
# dv = 36 # Number of design variables
# SB = np.hstack((lowBound,upBound)) # Stack upper and lower bounds


# lamHist, Lam, bestLam, Pi, Min, PAve, Ave = myGA(S,G,P,SB,lightboxSim,lightboxConstants)

# #%% Run best design to save data for plotting

# # Simulation constants
# lightboxConstants = {

#     # time-stepping constants
#     'stepFactor': 2.0, # factor by which to increase time step size
#     'timeMultiplier': 5.0, # multiplier to determine final time

#     # raytracing constants
#     'c': 3e8, # speed of light (m/s)
#     'rayDensity': 2000, # number of rays per face
#     'reflections': True, # marks if reflections are allowed
#     'rayTol': 0.000001, # threshold for ray absorption

#     'X1MinusOn': True, # marks if -x1 face will have rays
#     'X1PlusOn': True, # marks if +x1 face will have rays
#     'X2MinusOn': True, # marks if -x2 face will have rays
#     'X2PlusOn': True, # marks if +x2 face will have rays
#     'X3MinusOn': True, # marks if -x3 face will have rays
#     'X3PlusOn': True, # marks if +x3 face will have rays

#     # voxel constants
#     'deltaX1': 0.05, # voxel length (m)
#     'deltaX2': 0.05, # voxel width (m)
#     'deltaX3': 0.05, # voxel height (m)
#     'nVoxelX1': 20, # number of voxels in x1 direction
#     'nVoxelX2': 20, # number of voxels in x2 direction
#     'nVoxelX3': 20, # number of voxels in x3 direction

#     # target constants
#     'nTargetsX1': 9, # number of targets in x1 direction
#     'nTargetsX2': 3, # number of targets in x2 direction
#     'nTargetsX3': 3, # number of targets in x3 direction
#     'rackX1X3': True, # marks if rack is in x1-x3 plane
#     'rackX2X3': True, # marks if rack is in x2-x3 plane
#     'rackX1X2': True, # marks if rack is in x1-x2 plane
#     'rackX1': False, # marks if rack is in x1 direction
#     'rackX2': False, # marks if rack is in x2 direction
#     'rackX3': False, # marks if rack is in x3 direction

#     # power constants
#     'powerMagnitude': 1.0, # magnitude of power (to be multiplied by power dependency parameters)

#     # misc constants
#     'randseed': np.random.RandomState(201) # seed for ray variable randomizer

# }

# # run lightbox simulation
# Pi = lightboxSim(bestLam, lightboxConstants, movieSettings)
