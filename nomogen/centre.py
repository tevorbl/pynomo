#!/usr/bin/env python3

'''
verify that a good initial estimate is fornd when a
non vertical w axis passes thru (0.5,0.5)
'''

# pylint: disable=C

import sys
import math

import inspect
import os

sys.path.insert(0, "..")

from nomogen import Nomogen
from pynomo.nomographer import Nomographer

# get current file name
myfile = os.path.basename(inspect.stack()[0][1]).replace(".py", "")

# these allow combinations the v & w axes to be up or down
wup = vup = None


########################################
#
#  this is the target function,
#  - the function that the nomogram implements
#  - add the limits of the variables below
#
#  format is m = m(l,r), where l, m & r are respectively the values
#                        for the left, middle & right hand axes
#
#
########################################


# return value is the middle scale
# u & v are respectively the values on the left and right scales
#
# for this test, we require w(umin, vmax) == w(umax,vmin)
# so the lines intersect at the centre

# this is a sine curve that goes thru the centre
alpha = 0.25
scale = 1
def centre_sin(au,av):
    u = au  # math.sqrt(3*au + 1) - 1
    v = (av/scale)**2
    if not vup:
        v = 1-v
    f = lambda w: (u-v)*(0.5 - alpha*math.cos(math.pi*w)) + w - u
    df = lambda w: 1 + (u-v)*alpha*math.pi*math.sin(math.pi*w)
    w0 = (u+v)/2  # initial value
    n=0
    while abs(f(w0)) > 1e-13:
        w0 = w0 - f(w0)/df(w0)  # newton's method
        #print( f'{w0=}' )
        #print( f'{n=}' )
        if n > 20:
            sys.exit( 'centre_sin() could not find a value for w' )
        n += 1
    #print( f'{n=}' )
    if not wup:
        w0 = 1-w0
    return w0

# define the w axis as a parabola that touches the centre
# xw(w) = 0.5 + (w-0.5)^2
# so xw(0) = 0.5 + 1/4, xw(1) = 0.5 + 1/4,  xw(0.5) = 0.5
# determinant is
#  | 0  u  1 |
#  | xw w  1 |  = 0
#  | 1  v  1 |
#
# xw*(u-v) - (u-w) = 0
def centre_par(au,av):
    u = au
    v = av if vup else 1-av

    if math.isclose(u,v):
        yw = (u+v)/2      # +O((u-v)^2), also takes care of divide by zero
    else:
        d = math.sqrt(2*alpha*(u-v)*(u+v-1) + 1)
        yw = (alpha*(v-u) + 1 - d)/(2*alpha*(v-u))

   # check:
    xw = 0.5 + alpha*(yw-0.5)**2
    if not math.isclose(xw*(u-v), (u-yw), abs_tol=1e-5):
        print( f'par equation fail, u is {u}, v is {v}, w is {yw}')
        print( f'{xw*(u-v)=} , {(u-yw)=}' )

    return yw if wup else 1-yw



# range for the u scale (the left scale)
umin = 0
umax = 1

# range for the v scale (the right scale)
vmin = 0
vmax = scale


###############################################################
#
# nr points (Chebyshev nodes) needed to define the scales
# a higher value may be necessary if the scales are very non-linear
# a lower value is faster and makes a smoother curve,
#     but could be less accurate
NN = 4


##############################################
#
# definitions for the axes for pyNomo
# dictionary with key:value pairs

left_axis = {
    'u_min': umin,
    'u_max': umax,
    'title': 'u scale',
    'title_x_shift': 0.5,
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

right_axis = {
    'u_min': vmin,
    'u_max': vmax,
    'title': 'v scale',
    'title_x_shift': 0.5,
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

middle_axis = {
#    'u_min': wmin,
#    'u_max': wmax,
    'title': 'w scale',
    'title_x_shift': -0.2,
    'scale_type': 'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 2,
}

# assemble the above 3 axes into a block
block_params0 = {
    'block_type': 'type_9',
    'f1_params': left_axis,
    'f2_params': middle_axis,
    'f3_params': right_axis,

    # the isopleth connects the mid values of the outer axes
    # edit this for different values
    'isopleth_values': [[(left_axis['u_min'] + left_axis['u_max']) / 2, \
                         'x', \
                         (right_axis['u_min'] + right_axis['u_max']) / 2]],

    # log alignment errors
    # If this is missing or False then alignment error logs are disabled
    'LogAlignment': False,
}

# the nomogram parameters

main_params = {
#    'filename': myfile,
    'paper_height': 10,  # units are cm
    'paper_width': 10,
    'title_x': 7.0,
    'title_y': 2.0,
    'title_box_width': 8.0,
#    'title_str': r'$w = w(u,v)$',
    'block_params': [block_params0],

    # set the colour of the ispleth/index line here
    'isopleth_params': [ {'color': 'Red'}, ],

    'transformations': [('scale paper',)],

    'npoints': NN,

    # instead of forcing the ends of the axes to the corners of the unit square,
    # nomogen can shape the nomogram to minimise parallax errors
    # uncomment the following line to select this option
    #'muShape': 1,

    # text to appear at the foot of the nomogram
    # note that latex rules apply
    # a default string will appear if this is omitted
    # make this an empty string to have no footer text
    'footer_string': f'{myfile}: python {sys.version_info[0]}'
                                     f'.{sys.version_info[1]}'
                                     f'.{sys.version_info[2]}'
}

for func in [centre_par, centre_sin]:
    if func is centre_sin:
        fstr = 'sin'
    else:
        fstr = 'par'

    for vup in [True, False]:
        vstr = 'vup' if vup else 'vdown'
        for wup in [True, False]:
            wstr = 'wup' if wup else 'wdown'
            main_params['filename'] = myfile + f'.{vstr}.{wstr}.{fstr}.pdf'
            main_params['title_str'] = f'{vstr} {wstr} {fstr}'

            # range for the w scale (the middle scale)
            wvals = [func(uu,vv) for uu in [umin, umax] for vv in [vmin, vmax]]
            wmin = min(wvals)
            wmax = max(wvals)
            middle_axis['u_min'] = wmin
            middle_axis['u_max'] = wmax
            #print( f"wmin, wmax is {middle_axis['u_min']}, {middle_axis['u_max']}" )

            print( f"calculating the nomogram {vstr} {wstr} {fstr} ...")
            Nomogen(func, main_params)  # generate nomogram for the target function

            print("printing ", main_params['filename'], " ...")
            Nomographer(main_params)


            ### check where the w axis intercepts the top & bottom index lines

            if not math.isclose(wmin, 0, abs_tol = 1e-5):
                print( f'************ unexpected w min value, ({wmin}) ************ ' )
                continue

            if not math.isclose(wmax, 1, abs_tol = 1e-5):
                print( f'************ unexpected w max value, ({wmin}) ************ ' )
                continue

            if wup:
                t0 = middle_axis['f'](wmin) # the x coordinate of w at the bottom
                t1 = middle_axis['f'](wmax) # and at the top
            else:
                t0 = middle_axis['f'](wmax) # the x coordinate of w at the bottom
                t1 = middle_axis['f'](wmin) # and at the top

            if fstr == 'par':
                intercept0 = intercept1 = 0.5+alpha/4
            else: # 'sin'
                intercept0 = 0.5-alpha
                intercept1 = 0.5+alpha

            if not math.isclose( t0, intercept0, abs_tol=0.015 ):
                print( '*************** w axis intersects bottom index line,'
                       f' at {t0:.4g}, expected {intercept0:.4g} *************' )
            if not math.isclose( t1, intercept1, abs_tol=0.01 ):
                print( '*************** w axis intersects top index line,'
                       f' at {t1:.4g}, expected {intercept1:.4g} *************' )

            print()
