'''CD to the directory you want to run this in first'''

import os

path = os.getcwd()
filelist = os.listdir(path)
for fname in filelist:
    newfname = fname.replace('-',' - ')
    print fname, '  >>  ', newfname
    os.rename(fname, newfname)
