#!/usr/bin/env python3

'''
 nomogen test program
 verify orientation of axes, for all combinations of:
    v & w axis up/down
    linear/log scales
    symetrical/non-symetrical cases
    preScale on each axis

'''

# pylint: disable=C

import sys
import math

from pprint import  pprint

sys.path.insert(0, "..")

from nomogen import Nomogen
from pynomo.nomographer import Nomographer


# range for the u scale
umin = 1
umax = 5

# range for the v scale
vmin = 1
vmax = 5


def do_test( name_template ):
    ''' simple axes tests for each test case ... '''

    NN = 3        # default nr Chebyshev nodes

    for test_nr in range(1,17):
        test_nr += 1

        if test_nr == 1:
            w = lambda u,v: u+v
            scale_type = 'linear smart'
            title_str = r'$w = u + v$'
        elif test_nr == 2:
            w = lambda u, v: u - v
            scale_type = 'linear smart'
            title_str = r'$w = u - v$'
        elif test_nr == 3:
            w = lambda u, v: -(u + v)
            scale_type = 'linear smart'
            title_str = r'$w = -(u + v)$'
        elif test_nr == 4:
            w = lambda u, v: v - u
            scale_type = 'linear smart'
            title_str = r'$w = v - u$'
        elif test_nr == 5:
            w = lambda u, v: u*v
            scale_type = 'log smart'
            title_str = r'$w = uv$'
        elif test_nr == 6:
            w = lambda u, v: u / v
            scale_type = 'log smart'
            title_str = r'$w = {u \over v}$'
        elif test_nr == 7:
            w = lambda u, v: 1 / (u*v)
            scale_type = 'log smart'
            title_str = r'$w = {1 \over {uv}}$'
        elif test_nr == 8:
            w = lambda u, v: v / u
            scale_type = 'log smart'
            title_str = r'$w = {v \over u}$'
        elif test_nr == 9:
            w = lambda u,v: 0.4*u + v
            scale_type = 'linear smart'
            title_str = r'$w = 0.4u + v$'
        elif test_nr == 10:
            w = lambda u, v: 0.4*u - v
            scale_type = 'linear smart'
            title_str = r'$w = 0.4u-v$'
        elif test_nr == 11:
            w = lambda u, v: -(u + v*v)
            scale_type = 'linear smart'
            title_str = r'$w = -(u + v^2)$'
        elif test_nr == 12:
            w = lambda u, v: v - u*u
            scale_type = 'linear smart'
            title_str = r'$w = v - u^2$'
        elif test_nr == 13:
            w = lambda u, v: u**1.5*v
            scale_type = 'log smart'
            title_str = r'$w = u^{1.5}v$'
        elif test_nr == 14:
            w = lambda u, v: u**1.5 / v
            scale_type = 'log smart'
            title_str = r'$w = {u^{1.5} \over v}$'
        elif test_nr == 15:
            w = lambda u, v: 1 / (u*u*v)
            scale_type = 'log smart'
            title_str = r'$w = {1 \over {u^2v}}$'
        elif test_nr == 16:
            w = lambda u, v: v / u / u
            scale_type = 'log smart'
            title_str = r'$w = {v \over {u^2}}$'
        else:
            return #sys.exit( 'there is no test number ({})'.format(test_nr) )

        test_name = name_template.format(test_nr)


        t = test_nr - 1
        symStr = 'symmetrical' if t %16 < 8 else 'non-symmetrical'
        scaleStr = 'linear' if t % 8 < 4 else 'log'
        if name_template[-1] == 'a':
            wStr = 'down' if t % 4 < 2 else 'up'
            vStr = 'down' if t % 2 == 0 else 'up'
        else:
            wStr = 'up' if t % 4 < 2 else 'down'
            vStr = 'up' if t % 2 == 0 else 'down'
        tstr = '{}: {}, {}, w {}, v {}'.format(test_name, symStr, scaleStr, wStr, vStr)
        print( '\ntest', tstr )


        left_axis = {
            'u_min': umin,
            'u_max': umax,
            'title': r'$u$',
            'scale_type': scale_type,
            'tick_levels': 3,
            'tick_text_levels': 2,
        }
        if name_template[-1] == 'a':
            if 'log' in scale_type:
                left_axis['preScale'] = \
                    {'func': lambda t: math.log(wmin*wmax/t),
                     'inv': lambda t: wmin*wmax/math.exp(t)}
            else:
                left_axis['preScale'] = \
                    {'func': lambda t: (umin+umax)-t, 'inv': lambda t: (umin+umax)-t}
        elif 'preScale' in left_axis:
            left_axis.pop( 'preScale' )

        right_axis = {
            'u_min': vmin,
            'u_max': vmax,
            'title': r'$v$',
            'scale_type': scale_type,
            'tick_levels': 3,
            'tick_text_levels': 2,
        }
        if name_template[-1] == 'b':
            if 'log' in scale_type:
                right_axis['preScale'] = \
                    {'func': lambda t: math.log(wmin*wmax/t),
                     'inv': lambda t: wmin*wmax/math.exp(t)}
            else:
                right_axis['preScale'] = \
                    {'func': lambda t: (vmin+vmax-t),
                     'inv': lambda t: (vmin+vmax)-t}
        elif 'preScale' in right_axis:
            right_axis.pop( 'preScale' )

        # automagically get the w scale range
        wtmp = [ w(umin, vmin), w(umax, vmin), w(umin, vmax), w(umax, vmax) ]
        wmin = min(wtmp)
        wmax = max(wtmp)
        middle_axis = {
            'u_min': wmin,
            'u_max': wmax,
            'title': r'$w$',
            'scale_type': scale_type,
            'tick_levels': 3,
            'tick_text_levels': 2,
        }
        if name_template[-1] == 'c':
            if 'log' in scale_type:
                middle_axis['preScale'] = \
                    {'func': lambda t: math.log(wmin*wmax/t),
                     'inv': lambda t: wmin*wmax/math.exp(t)}
            else:
                middle_axis['preScale'] = {'func': lambda t: (wmin+wmax)-t,
                                               'inv': lambda t: (wmin+wmax)-t}
        elif 'preScale' in middle_axis:
            middle_axis.pop( 'preScale' )

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
        }

        # the nomogram parameters
        main_params = {
            'filename': test_name,
            'paper_height': 10,  # units are cm
            'paper_width': 10,
            'title_x': 5.0 if test_nr in [9, 10, 11, 12] else 7.0,
            'title_y': 2.0,
            'title_box_width': 8.0,
            'title_str': title_str,
            'block_params': [block_params0],
            'transformations': [('scale paper',)],
            'isopleth_params': [{'color': 'Red'}],
            'npoints': NN,

            # text to appear at the foot of the nomogram
            # make this null string for nothing
            # a default string will appear if this is omitted
            'footer_string': tstr,
            'trace': 'trace_init'
        }



        #pprint( main_params )

        print("calculating the nomogram ...")
        try:
            Nomogen(w, main_params)  # generate nomogram for the target function
            main_params['filename'] += '.pdf'
            print("printing ", main_params['filename'], " ...")
            Nomographer(main_params)
        except ValueError:
            print( test_name, 'failed')

    # end of do_test()


do_test( 'axestest.{}' )  # normal
do_test( 'axestest.{}a' )  # preScale on left axis
do_test( 'axestest.{}b' )  # preScale on right axis
do_test( 'axestest.{}c' )  # preScale on middle axis

################# end of axestest.py ###############


