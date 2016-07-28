#  coding = utf-8

import os
 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.


home_dir = os.getenv("HOME")
working_dir = str(home_dir) + '/hcss/workspace/POC/'
pool_dir = str(working_dir) + 'pools/'
list_camera = ['red', 'blue']

obsids = {}
obsids[1342208884] = "SEDA"
obsids[1342208885] = "SEDB"
obsids[1342244919] = "SEDC"
obsids[1342188037] = "SEDD"
obsids[1342188038] = "SEDE"

nLowFreq = 4
smoothing = 'filter'
gaussianFilterWidth = 50
medianFilterWidth = 15

#
# First observation
#
# obsid 1342208884 - camera red - normalization
camera = 'red'
sey = '1342208884_red_norm'
obs_1342208884_R_N = getObservation(obsid=1342208884,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342208884_R_N)    

slicedFinalCubes = obs_1342208884_R_N.level2.red.rcube.product
#for slice in range(len(slicedFinalCubes.refs)):
print len(slicedFinalCubes.refs)
c1_1342208884_rnc, c9_1342208884_rnc, c129_1342208884_rnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208884 - camera blue - normalization
camera = 'blue'
sey = '1342208884_blue_norm'
obs_1342208884_B_N = getObservation(obsid=1342208884,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342208884_B_N)  

slicedFinalCubes = obs_1342208884_B_N.level2.blue.rcube.product
print len(slicedFinalCubes.refs)
c1_1342208884_bnc, c9_1342208884_bnc, c129_1342208884_bnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208884 - camera red - poc 5
camera = 'red'
sey = '1342208884_red_corr_5'
obs_1342208884_R_POC_5 = getObservation(obsid=1342208884,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208884_R_POC_5)    

slicedFinalCubes = obs_1342208884_R_POC_5.level2.red.rcube.product
print len(slicedFinalCubes.refs)
c1_1342208884_rpoc5, c9_1342208884_rpoc5, c129_1342208884_rpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208884 - camera blue - poc 5
camera = 'blue'
sey = '1342208884_blue_corr_5'
obs_1342208884_B_POC_5 = getObservation(obsid=1342208884,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208884_B_POC_5)    

slicedFinalCubes = obs_1342208884_B_POC_5.level2.blue.rcube.product
c1_1342208884_bpoc5, c9_1342208884_bpoc5, c129_1342208884_bpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208884 - camera red - poc 9
camera = 'red'
sey = '1342208884_red_corr_9'
obs_1342208884_R_POC_9 = getObservation(obsid=1342208884,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208884_R_POC_9)    

slicedFinalCubes = obs_1342208884_R_POC_9.level2.red.rcube.product
c1_1342208884_rpoc9, c9_1342208884_rpoc9, c129_1342208884_rpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208884 - camera blue - poc 9
camera = 'blue'
sey = '1342208884_blue_corr_9'
obs_1342208884_B_POC_9 = getObservation(obsid=1342208884,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208884_B_POC_9)    

slicedFinalCubes = obs_1342208884_B_POC_9.level2.blue.rcube.product
c1_1342208884_bpoc9, c9_1342208884_bpoc9, c129_1342208884_bpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)
"""
#
# Second observation
#
# obsid 1342208885 - camera red - normalization
camera = 'red'
sey = '1342208885_red_norm'
obs_1342208885_R_N = getObservation(obsid=1342208885,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342208885_R_N)    

slicedFinalCubes = obs_1342208885_R_N.level2.red.rcube.product
c1_1342208885_rnc, c9_1342208885_rnc, c129_1342208885_rnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208885 - camera blue - normalization
camera = 'blue'
sey = '1342208885_blue_norm'
obs_1342208885_B_N = getObservation(obsid=1342208885,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342208885_B_N)  

slicedFinalCubes = obs_1342208885_B_N.level2.blue.rcube.product
c1_1342208885_bnc, c9_1342208885_bnc, c129_1342208885_bnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208885 - camera red - poc 5
camera = 'red'
sey = '1342208885_red_corr_5'
obs_1342208885_R_POC_5 = getObservation(obsid=1342208885,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208885_R_POC_5)    

slicedFinalCubes = obs_1342208885_R_POC_5.level2.red.rcube.product
c1_1342208885_rpoc5, c9_1342208885_rpoc5, c129_1342208885_rpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208885 - camera blue - poc 5
camera = 'blue'
sey = '1342208885_blue_corr_5'
obs_1342208885_B_POC_5 = getObservation(obsid=1342208885,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208885_B_POC_5)    

slicedFinalCubes = obs_1342208885_B_POC_5.level2.blue.rcube.product
c1_1342208885_bpoc5, c9_1342208885_bpoc5, c129_1342208885_bpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208885 - camera red - poc 9
camera = 'red'
sey = '1342208885_red_corr_9'
obs_1342208885_R_POC_9 = getObservation(obsid=1342208885,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208885_R_POC_9)    

slicedFinalCubes = obs_1342208885_R_POC_9.level2.red.rcube.product
c1_1342208885_rpoc9, c9_1342208885_rpoc9, c129_1342208885_rpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342208885 - camera blue - poc 9
camera = 'blue'
sey = '1342208885_blue_corr_9'
obs_1342208885_B_POC_9 = getObservation(obsid=1342208885,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342208885_B_POC_9)    

slicedFinalCubes = obs_1342208885_B_POC_9.level2.blue.rcube.product
c1_1342208885_bpoc9, c9_1342208885_bpoc9, c129_1342208885_bpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

#
# Third observation
#
# obsid 1342244919 - camera red - normalization
camera = 'red'
sey = '1342244919_red_norm'
obs_1342244919_R_N = getObservation(obsid=1342244919,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342244919_R_N)    

slicedFinalCubes = obs_1342244919_R_N.level2.red.rcube.product
c1_rnc, c9_rnc, c129_rnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                  smoothing=smoothing,
                                                  width=gaussianFilterWidth,
                                                  preFilterWidth=medianFilterWidth,
                                                  nLowFreq=nLowFreq,
                                                  calTree=calTree,
                                                  verbose=False)

# obsid 1342244919 - camera blue - normalization
camera = 'blue'
sey = '1342244919_blue_norm'
obs_1342244919_B_N = getObservation(obsid=1342244919,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342244919_B_N)  

slicedFinalCubes = obs_1342244919_B_N.level2.blue.rcube.product
c1_1342244919_bnc, c9_1342244919_bnc, c129_1342244919_bnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342244919 - camera red - poc 5
camera = 'red'
sey = '1342244919_red_corr_5'
obs_1342244919_R_POC_5 = getObservation(obsid=1342244919,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342244919_R_POC_5)    

slicedFinalCubes = obs_1342244919_R_POC_5.level2.red.rcube.product
c1_1342244919_rpoc5, c9_1342244919_rpoc5, c129_1342244919_rpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342244919 - camera blue - poc 5
camera = 'blue'
sey = '1342244919_blue_corr_5'
obs_1342244919_B_POC_5 = getObservation(obsid=1342244919,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342244919_B_POC_5)    

slicedFinalCubes = obs_1342244919_B_POC_5.level2.blue.rcube.product
c1_1342244919_bpoc5, c9_1342244919_bpoc5, c129_1342244919_bpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342244919 - camera red - poc 9
camera = 'red'
sey = '1342244919_red_corr_9'
obs_1342244919_R_POC_9 = getObservation(obsid=1342244919,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342244919_R_POC_9)    

slicedFinalCubes = obs_1342244919_R_POC_9.level2.red.rcube.product
c1_1342244919_rpoc9, c9_1342244919_rpoc9, c129_1342244919_rpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342244919 - camera blue - poc 9
camera = 'blue'
sey = '1342244919_blue_corr_9'
obs_1342244919_B_POC_9 = getObservation(obsid=1342244919,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342244919_B_POC_9)    

slicedFinalCubes = obs_1342244919_B_POC_9.level2.blue.rcube.product
c1_1342244919_bpoc9, c9_1342244919_bpoc9, c129_1342244919_bpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)
"""
#
# Forth observation
#
# obsid 1342188037 - camera red - normalization
camera = 'red'
sey = '1342188037_red_norm'
obs_1342188037_R_N = getObservation(obsid=1342188037,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342188037_R_N)    

slicedFinalCubes = obs_1342188037_R_N.level2.red.rcube.product
c1_1342188037_rnc, c9_1342188037_rnc, c129_1342188037_rnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208884 - camera blue - normalization
camera = 'blue'
sey = '1342188037_blue_norm'
obs_1342188037_B_N = getObservation(obsid=1342188037,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342188037_B_N)  

slicedFinalCubes = obs_1342188037_B_N.level2.blue.rcube.product
c1_1342188037_bnc, c9_1342188037_bnc, c129_1342188037_bnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342208884 - camera red - poc 5
camera = 'red'
sey = '1342188037_red_corr_5'
obs_1342188037_R_POC_5 = getObservation(obsid=1342188037,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188037_R_POC_5)    

slicedFinalCubes = obs_1342188037_R_POC_5.level2.red.rcube.product
c1_1342188037_rpoc5, c9_1342188037_rpoc5, c129_1342188037_rpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188037 - camera blue - poc 5
camera = 'blue'
sey = '1342188037_blue_corr_5'
obs_1342188037_B_POC_5 = getObservation(obsid=1342188037,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188037_B_POC_5)    

slicedFinalCubes = obs_1342188037_B_POC_5.level2.blue.rcube.product
c1_1342188037_bpoc5, c9_1342188037_bpoc5, c129_1342188037_bpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188037 - camera red - poc 9
camera = 'red'
sey = '1342188037_red_corr_9'
obs_1342188037_R_POC_9 = getObservation(obsid=1342188037,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188037_R_POC_9)    

slicedFinalCubes = obs_1342188037_R_POC_9.level2.red.rcube.product
c1_1342188037_rpoc9, c9_1342188037_rpoc9, c129_1342188037_rpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188037 - camera blue - poc 9
camera = 'blue'
sey = '1342188037_blue_corr_9'
obs_1342188037_B_POC_9 = getObservation(obsid=1342188037,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188037_B_POC_9)    

slicedFinalCubes = obs_1342188037_B_POC_9.level2.blue.rcube.product
c1_1342188037_bpoc9, c9_1342188037_bpoc9, c129_1342188037_bpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

#
# Fifth observation
#
# obsid 1342188038 - camera red - normalization
camera = 'red'
sey = '1342188038_red_norm'
obs_1342188038_R_N = getObservation(obsid=1342188038,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342188038_R_N)    

slicedFinalCubes = obs_1342188038_R_N.level2.red.rcube.product
c1_1342188038_rnc, c9_1342188038_rnc, c129_1342188038_rnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342188038 - camera blue - normalization
camera = 'blue'
sey = '1342188038_blue_norm'
obs_1342188038_B_N = getObservation(obsid=1342188038,
                                    poolLocation=str(pool_dir), poolName=sey)
calTree = getCalTree(obs=obs_1342188038_B_N)  

slicedFinalCubes = obs_1342188038_B_N.level2.blue.rcube.product
c1_1342188038_bnc, c9_1342188038_bnc, c129_1342188038_bnc = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                   smoothing=smoothing,
                                                                                   width=gaussianFilterWidth,
                                                                                   preFilterWidth=medianFilterWidth,
                                                                                   nLowFreq=nLowFreq,
                                                                                   calTree=calTree,
                                                                                   verbose=False)

# obsid 1342188038 - camera red - poc 5
camera = 'red'
sey = '1342188038_red_corr_5'
obs_1342188038_R_POC_5 = getObservation(obsid=1342188038,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188038_R_POC_5)    

slicedFinalCubes = obs_1342188038_R_POC_5.level2.red.rcube.product
c1_1342188038_rpoc5, c9_1342188038_rpoc5, c129_1342188038_rpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188038 - camera blue - poc 5
camera = 'blue'
sey = '1342188038_blue_corr_5'
obs_1342188038_B_POC_5 = getObservation(obsid=1342188038,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188038_B_POC_5)    

slicedFinalCubes = obs_1342188038_B_POC_5.level2.blue.rcube.product
c1_1342188038_bpoc5, c9_1342188038_bpoc5, c129_1342188038_bpoc5 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188038 - camera red - poc 9
camera = 'red'
sey = '1342188038_red_corr_9'
obs_1342188038_R_POC_9 = getObservation(obsid=1342188038,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188038_R_POC_9)    

slicedFinalCubes = obs_1342188038_R_POC_9.level2.red.rcube.product
c1_1342188038_rpoc9, c9_1342188038_rpoc9, c129_1342188038_rpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)

# obsid 1342188038 - camera blue - poc 9
camera = 'blue'
sey = '1342188038_blue_corr_9'
obs_1342188038_B_POC_9 = getObservation(obsid=1342188038,
                                        poolLocation=str(pool_dir),
                                        poolName=sey)
calTree = getCalTree(obs=obs_1342188038_B_POC_9)    

slicedFinalCubes = obs_1342188038_B_POC_9.level2.blue.rcube.product
c1_1342188038_bpoc9, c9_1342188038_bpoc9, c129_1342188038_bpoc9 = extractCentralSpectrum(slicedFinalCubes.get(0),
                                                                                         smoothing=smoothing,
                                                                                         width=gaussianFilterWidth,
                                                                                         preFilterWidth=medianFilterWidth,
                                                                                         nLowFreq=nLowFreq,
                                                                                         calTree=calTree,
                                                                                         verbose=False)


plot_1342208884_c9 = openSE(c9_1342208884_rnc, display=1)
plot_1342208884_c9.add(c9_1342208884_rpoc5)
plot_1342208884_c9.add(c9_1342208884_rpoc9)
