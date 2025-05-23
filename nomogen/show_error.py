#!/usr/bin/env python3

'''
simple program to plot nomogen alignment errors using matplotlib
argument, which is the filename of the alignment errors creatyed by nomogen
'''

# relevant only for python >= 3.6
# pylint: disable=consider-using-f-string

import sys
import os
import shutil

import matplotlib.pyplot as plt

import numpy as np

if len(sys.argv) < 2:
    print( 'an alignment error file is needed' )
    print(
        '''\nUsage:
            {} filename
            where filename is name of alignment error generated by nomogen'''.format(sys.argv[0])
    )
    sys.exit('quitting')


fn = sys.argv[1]
if os.path.exists(fn):
    shutil.copyfile(fn, 'nomo_error.py' )
else:
    sys.exit("The alignment error file '{}' does not exist - quitting".format(fn))

import nomo_error

# no error checking
print( 'using nomo_error file', nomo_error.nomo_id )
ax = np.array([x['w'] for x in nomo_error.Nomo_error])
ay = np.array([y['error'] for y in nomo_error.Nomo_error])
if np.max(ay) >= 0.1 or np.min(ay) <= -0.1:
    scale = 'mm'
else:
    ay *= 1000
    scale = 'microns'

# the plot
efig, eplt = plt.subplots()

eplt.scatter(ax, ay, color='red', s=5, marker='.' )

eplt.grid(True, linestyle='dashed', )

# embolden the zero error line
y_ticks = eplt.get_yticks(  )
ind = np.where( y_ticks==0)[0][0]  # find index where error == 0 mm

gridlines = eplt.yaxis.get_gridlines()
gridlines[ind].set_color("black")
gridlines[ind].set_linewidth(1.5)

plt.title( 'Alignment Error for {} nomogram'.format(nomo_error.nomo_id['file']) )
plt.xlabel( '{} scale'.format(nomo_error.nomo_id['title']) )
plt.ylabel('error {}'.format(scale))

try:
    plt.savefig( fn.replace('py', 'pdf') )
except ValueError:
    # maybe matplotlib is unable to parse the title of the w scale
    plt.xlabel( 'w scale' )  # use backup title
    plt.savefig( fn.replace('py', 'pdf') )

plt.show()
