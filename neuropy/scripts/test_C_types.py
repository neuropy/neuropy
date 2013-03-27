"""Test C_types.pyx"""

import pyximport
pyximport.install(build_in_temp=False, inplace=True)
import C_types # .pyx file

C_types.test()
