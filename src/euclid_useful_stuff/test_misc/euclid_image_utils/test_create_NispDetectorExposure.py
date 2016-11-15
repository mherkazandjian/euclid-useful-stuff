"""
<keywords>
test, python, euclid, image, utils, nisp, detector, exposure, imagef
</keywords>
<description>
example for creating an 3 layer nisp detector exposure object
</description>
<seealso>
</seealso>
"""
from __future__ import print_function

from numpy import zeros

from EuclidImageBinding import NispDetectorExposure, ImageF, MaskUI

im_size = 2040

# define the arrays of each layer
sci_array = ImageF(zeros((im_size, im_size), 'f4'))
var_array = ImageF(zeros((im_size, im_size), 'f4'))
dq_array = MaskUI(zeros((im_size, im_size), 'uint32'))

# create the detector object from the 3 layers
det_exposure = NispDetectorExposure(sci_array, var_array, dq_array)

print('done')
