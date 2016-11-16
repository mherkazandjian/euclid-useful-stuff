"""
<keywords>
test, python, euclid, image, utils, exposure, mef
</keywords>
<description>
example of read the survey exposure file (three layers)
</description>
<seealso>
</seealso>
"""
from __future__ import print_function

import os
from os.path import join

from EuclidImageBinding import NispSurveyExposure

fname = os.path.expanduser('~/euclid/data/euclid/simulated/SC2/NIP_R3/NIP_R3_3L/EUC-TEST-Y-2016-04-26T172918.532Z.fits')

print('input fits file:\n\t%s' % fname)

exposure = NispSurveyExposure.readFits(fname)

detector = exposure.getDetector(0)

metadata = detector.getMetadata()

gain = metadata.get('GAIN')
readnoise = metadata.get('RDNOISE')
saturate = metadata.get('SATURATE')

sci = detector.getScience().getArray()
print('mean pixel value for the science layer  = ', sci.mean())

var = detector.getVariance().getArray()
print('mean pixel value for the variance layer = ', var.mean())

mask = detector.getMask().getArray()
print('mean pixel value for the variance layer = ', var.mean())

print('done')
