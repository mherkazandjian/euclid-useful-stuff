"""
<keywords>
test, python, euclid, image, utils, property, list
</keywords>
<description>
example for creating a property list
</description>
<seealso>
</seealso>
"""
from __future__ import print_function

from numpy import zeros

from EuclidImageBinding import PropertyList

plist = PropertyList()
plist.add('xxx', -666)
print(plist)

print('done')
