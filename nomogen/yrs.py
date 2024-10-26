#!/usr/bin/env python3

#nomogen example program

import sys

import inspect
import os

sys.path.insert(0, "..")

from math import *

from nomogen import Nomogen
from pynomo.nomographer import Nomographer

# get current file name
myfile = os.path.basename(inspect.stack()[0][1]).replace(".py", "")

# alternative with no external dependencies - it works most of the time
#  myfile =  __name__ == "__main__" and (__file__.endswith(".py") and __file__.replace(".py", "") or "nomogen")
#             or __name__,


########################################
#
#  this is the target function,
#  - the limits of the variables
#  - the function that the nonogram implements
#
#  format is m = m(l,r), where l, m & r are respectively the values
#                        for the left, middle & right hand axes
########################################

# nr years to achieve future value fv by investing $1 every year
# - years is return value and will be the middle scale,
# - first argument is fv and will become the left scale, and
# - second argument is r and will become the right scale
#
#   fv = ((1 + i)**y - 1) / i
#
def yrs(fv, r):
    i = r/100
    y = log(fv*i + 1)/log(i+1)
    if abs((((1 + i) ** y) - 1)/i - fv) > 1e-10:
        print("yrs: equation fault, y = ", y, ", i is ", i, ", fv is ", fv)
        sys.exit("quitting")
    return y

fvmin = 1;   fvmax = 33     # required future value
imin = 0.5; imax = 5        # interest rates
ymin = yrs(fvmin, imax)    # nr years
ymax = yrs(fvmax, imin)

###############################################################
#
# nr Chebyshev nodes needed to define the scales
# a higher value may be necessary if the scales are very non-linear
# a lower value increases speed, makes a smoother curve, but could introduce errors
NN = 5


##############################################
#
# definitions for the axes for pyNomo
# dictionary with key:value pairs

left_axis = {
    'u_min': fvmin,
    'u_max': fvmax,
    'title': r'future value',
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

right_axis = {
    'u_min': imin,
    'u_max': imax,
    'title': r'interest rate',
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

middle_axis = {
    'u_min': ymin,
    'u_max': ymax,
    'title': r'years',
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

block_params0 = {
    'block_type': 'type_9',
    'f1_params': left_axis,
    'f2_params': middle_axis,
    'f3_params': right_axis,
    'transform_ini': False,
    'isopleth_values': [[(left_axis['u_min'] + left_axis['u_max'])/2, \
                         'x', \
                         (right_axis['u_min'] + right_axis['u_max'])/2]],

    # log alignment errors
    # If this is missing or False then alignment error logs are disabled
    'LogAlignment': True,
}

main_params = {
    'filename': myfile,
    'paper_height': 10, # units are cm
    'paper_width': 10,
    'title_x': 6.0,
    'title_y': 1.0,
    'title_box_width': 8.0,
    'title_str':r'future value of \$1 invested each year',
    'extra_texts':[
        {'x':5,
         'y':2,
         'text':r'$FV = {(1+i)^y-1 \over i}$',
         'width':5,
         }],
    'block_params': [block_params0],
    'transformations': [('scale paper',)],
    'npoints': NN
}

print("calculating the nomogram ...")
Nomogen(yrs, main_params)  # generate nomogram for yrs function

main_params['filename'] += '.pdf'
print("printing ", main_params['filename'], " ...")
Nomographer(main_params)
