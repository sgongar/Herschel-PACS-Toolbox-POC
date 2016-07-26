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

obsids = {}
obsids[1342208884] = "SEDA"
obsids[1342208885] = "SEDB"
obsids[1342244919] = "SEDC"
obsids[1342188037] = "SEDD"
obsids[1342188038] = "SEDE"

for i in range(len(obsids.keys())):
    obs = getObservation(obsids.keys()[i], useHsa = 1)

    execfile(str(working_dir) + 'ChopNodBackgroundNormalizationRange.py')
    name = str(obsids.key()[i]) + 'norm'
    saveObservation(obs, poolLocation = pool_dir, poolName = name)

    execfile(str(working_dir) + 'ChopNodRangeScanPointingCorrection_5.py')
    name = str(obsids.keys()[i]) + 'corr_5'
    saveObservation(obs, poolLocation = pool_dir, poolName = name)

    execfile(str(working_dir) + 'ChopNodRangeScanPointingCorrection_9.py')
    name = str(obsids.keys()[i]) + 'corr_9'
    saveObservation(obs, poolLocation = pool_dir, poolName = name)