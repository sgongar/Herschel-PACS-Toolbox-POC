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

verbose = False

obsids = {}
# obsids[1342208884] = "SEDA"
obsids[1342208885] = "SEDB"
obsids[1342244919] = "SEDC"
# obsids[1342188037] = "SEDD"
# obsids[1342188038] = "SEDE"

list_camera = ['red', 'blue']

def plot_graph(first, second, three):
    print first
    print second
    print three

for j in range(len(obsids.keys())):
    for k in range(len(list_camera)):
        camera = list_camera[k]
        
        obs = getObservation(obsids.keys()[j], useHsa=1)
        """
        print "ALLES GUT ZUM 1"
        execfile(str(working_dir) + 'ChopNodBackgroundNormalizationRange.py')
        name = str(obsids.keys()[j]) + '_' + str(camera) + '_norm'
        saveObservation(obs, poolLocation = pool_dir, poolName=name)

        print "ALLES GUT ZUM 2"
        execfile(str(working_dir) + 'ChopNodRangeScanPointingCorrection_5.py')
        name = str(obsids.keys()[j]) + '_' + str(camera) + '_corr_5'
        saveObservation(obs, poolLocation=pool_dir, poolName=name)
        """
        print "ALLES GUT ZUM 3"
        execfile(str(working_dir) + 'ChopNodRangeScanPointingCorrection_9.py')
        name = str(obsids.keys()[j]) + '_' + str(camera) + '_corr_9'
        saveObservation(obs, poolLocation=pool_dir, poolName=name)
"""
# Produce a plot after each iteration
for j in range(len(obsids.keys())):
    for k in range(len(list_camera)):
        camera = list_camera[k]
 
        sey = str(obsids.keys()[j]) + '_' + str(camera) + '_norm'
        obs_case1 = getObservation(obsid=obsids.keys()[j],
                                   poolLocation=str(pool_dir), poolName=sey)
        sey = str(obsids.keys()[j]) + '_' + str(camera) + '_corr_5'
        obs_case2 = getObservation(obsid=obsids.keys()[j],
                                   poolLocation=str(pool_dir), poolName=sey)

        sey = str(obsids.keys()[j]) + '_' + str(camera) + '_corr_9'
        obs_case3 = getObservation(obsid=obsids.keys()[j],
                                   poolLocation=str(pool_dir), poolName=sey)
        
        # plot_graph('blabla1', 'blabla2', 'blabla3')\
"""