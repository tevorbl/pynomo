#!/usr/bin/env python3

"""
    nomogen.py

    auto generate a 3 line nomogram

    given a function w=w(u,v), numerically derive a nomogram
    generate the pdf with pynomo


    Copyright (C) 2021-2026  Trevor Blight

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""



import sys
import datetime
import functools

import math
import numpy as np
import scipy.interpolate


# relevant only for python >= 3.6
# pylint: disable=consider-using-f-string

# pylint: disable=invalid-name
# pylint: max-line-length=128


###############################################
#
# global variables
# used by by optimiser to manage iterations and report progress

itNr = 0                           # iteration counter
eAcc = eDeru = eDerv = eShape = 0  # total alignment error, derivative errors & shape error
old_cost = cur_cost = 0            # current cost & cost of previous iteration
lcost1 = 0                          # log cost at first iteration


def Nomogen( func, main_params ):
    """
    Nomogen - numerically determine a nomogram from an equation
    inputs:
             func(a,b), return value is middle scale,
                        a & b are values on left and right scales
             params     as in Nomographer, but with these additions

                        'npoints' :
                         a measure of how many points are needed to
                         define each scale line in the nomogram

                        'muShape' :
                        set to non-zero so the ends of the outer axes
                        are not fixed to the corners of the unit square

                        'LogAlignment':
                        if set True write a file containing alignment error
                        measurements.

                        'footer_string':
                        text to appear at the foot of the nomogram.
                        Default text appears if the is missing.
    """


    ###################################
    #
    # first check for valid inputs
    #
    ###################################

    if main_params['block_params'][0]['block_type'] != 'type_9':
        sys.exit( '''
        'type_9' block expected, found '{}'
        '''.format(main_params['block_params'][0]['block_type']) )

    # NN = nr Chebyshev nodes for the scales
    # the nodes have an index in the range 0 .. NN-1
    if 'pdegree' in main_params:
        pstr = 'pdegree'
        NN = main_params['pdegree']
        print( '''
        WARNING: parameter 'pdegree' has been renamed 'npoints' to more accurately describe its use.
        Please update your code.
        ''' )
    else:
        pstr = 'npoints'
        if 'npoints' in main_params:
            NN = main_params['npoints']
        else:
            print( "\n'npoints' parameter not provided,", end=' ' )
            NN = 9

    if not isinstance(NN, int):
        print( "\nparameter '{}' must be an integer, found {} - ".format(pstr, type(NN)), end=' ' )
        NN = 9
    else:
        if NN < 3:
            sys.exit( "{} must be >= 3 - quitting".format(pstr) )
    print( "using", NN, "Chebyshev points" )

    if 'muShape' in main_params and main_params['muShape']:
        # magic scaling constant
        # tune this so shape error eShape is same order as
        # or slightly smaller than eAcc & other costs
        muShape = 1.0e-9
        print( 'setting nomogram for best measureability, muShape is', muShape )
    else:
        muShape = 0


    width = 10 * main_params['paper_width']  # convert cm -> mm
    height = 10 * main_params['paper_height']

    # allow nomogram tolerance == +/- 0.1 mm
    # this is how many dots 0.1mm square fit into the nomogram area
    resolution = 10 * 10 * width * height

    # grid for NN chebyshev nodes scaled to the interval [0,1]
    # nodes must be in increasing order
    chebGrid = (1 - np.cos(np.linspace(0, math.pi, NN)))/2
    # print( chebGrid )

    # enable logging
    # undocumented feature, for developers only
    # "trace_init", "trace_cost", "trace_alignment", "trace_result", "trace_derivative"
    if 'trace' in main_params:
        trace = main_params['trace']
        print('trace is ', trace)
    else:
        trace = {  }


    #############################################################
    #
    sqrt_eps = math.sqrt( sys.float_info.epsilon )
    def derivative(f, x, x1, x0):
        '''
        find the derivative of function f at position x
        detect range and value errors
        x1 & x0 are max & min values of x, used to calculate step size

        formulas taken from Fornberg:
        "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
        Mathematics of Computation, Vol 51, nr 184, Oct 1988, pp 699-706
        '''

        # scale by 1024 for O(h^4)
        # higher order means larger step size is OK, leaving more bits for accuracy
        # use O(h^4) methods, seems to give 7 or 8 digits of accuracy in practice
        step_size = sqrt_eps * 1024
        h = (x1-x0)*step_size   # step size

        try:
            r = (f(x-2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h)  # O(h^4)

        except ValueError as e:
            # try one sided O(h^4) calculation
            if x >= x1-2*h:
                r = (-25*f(x) + 48*f(x-h) - 36*f(x-2*h) + 16*f(x-3*h) - 3*f(x-4*h))/(12*h)
                if 'trace_derivative' in trace:
                    print( f'attempting to evaluate above valid range ({x1}), r is {r}' )
            elif x <= x0+2*h:
                r = (-25*f(x) + 48*f(x+h) - 36*f(x+2*h) + 16*f(x+3*h) - 3*f(x+4*h))/(12*h)
                if 'trace_derivative' in trace:
                    print( f'attempting to evaluate below valid range ({x0}), r is {r}' )
            else:
                # cannot recover, giving up
                r =  math.nan
                print( '-----------------------------------' )
                print( f'derivative error ({e}) at x={x}, max={x1}, min={x0}' )
        return r


    class Axis:
        '''
        encapsulate an axis of the nomogram
        axis equation, x & y coordinates
        Chebyshev grid
        preScale function
        '''

        # constructor
        def __init__( self, aparams, limits=None ):
            '''
            axis constructor
            aparams is the axis parameter block
            limits is a list of values that may limit the range of the axis
               this might be used for a w axis that has limits determined by
               function evaluations
            get the min & max input values
            check for sanity
            plotted values are linear, log or user defined by 'preScale' parameter
            convert between plotted and user values
            if log scale, use log of corresponding variable
            define plot * user conversion functions
            set up the Chebyshev nodes and differentiation matrix
            '''
            self.params = aparams
            self.title = aparams['title']
            if limits is None:
                self.max_user = aparams['u_max']
                self.min_user = aparams['u_min']
                if self.max_user < self.min_user:
                    sys.exit("error in {} scale: umax ({}) is less than umin ({})"
                             .format(self.title, self.max_user, self.min_user))
            else:
                a = limits + [aparams['u_max'], aparams['u_min']]
                self.max_user = max(a)
                self.min_user = min(a)

            if 'preScale' in aparams:
                self.plot = aparams['preScale']['func']
                self.user = aparams['preScale']['inv']
                t = (self.min_user + self.max_user)/2
                if not math.isclose( t, self.user(self.plot(t)) ):
                    sys.exit( f"""
preScale functions 'func' and 'inv' are not inverses of each other
func({t}) is {self.plot(t)}, inv({self.plot(t)}) is {self.user(t)}
""" )
            elif 'scale_type' in aparams and 'log' in aparams['scale_type']:
                if self.min_user <= 0:
                    sys.exit("error in {} scale: cannot have min ({}) <= 0 for log scale".
                             format(self.title, self.min_user))
                self.plot = math.log
                self.user = math.exp
            else:
                self.user = self.plot = lambda t: t

            # convert from user to plot values
            a = self.plot(self.min_user)
            b = self.plot(self.max_user)
            self.min = min(a,b)
            self.max = max(a,b)

            # with min & max finalised, we can add nodes
            # scale to the axis
            # these arrays define the x & y coordinates of the scales
            self.nodes = self.min + (self.max - self.min) * chebGrid
            self.bary = scipy.interpolate.BarycentricInterpolator(self.nodes)
            self.diffmat = self._cheb()  # the differentiation matrix


        def _cheb(self):
            '''
            compute the derivative matrix for a chebyshev grid
            see cheb.m in "Spectral Methods in Matlab", Ch 6, by Lloyd N. Trefethen
            x: the chebyshev grid
            return the derivative matrix
            '''

            n = len(self.nodes)  # nomogen guarantees that n >= 3

            # replicate the Chebyshev grid onto a nxn matrix
            XX = np.tile( self.nodes, (n,1) )

            # coefficients, alternating 1/-1/..., first & last entries are +/-2
            c = np.empty((n,1))  # column vector
            c[::2] = 1            # evens
            c[1::2] = -1          # odds
            c[0] = 2; c[-1] *= 2  # first & last

            D = c * (1/c).T  # nxn matrix

            # off diagonal entries
            dXX = XX.T - XX + np.eye(n)
            D = D / dXX

            # diagonal entries - subtract sum of each row
            for ii in range(n):
                D[ii][ii] -= D[ii].sum()

            return D


        def __call__( self, a, afn ):
            '''
            evaluate function at a using barycentric interpolation,
            given the objects array of Chebyshev nodes and function values
            a   - perform barycentric interpolation at position a
            afn - array of function values at this objects nodes
            '''

            # special case: a is one of the grid points
            if a in self.nodes:
                t = np.searchsorted(self.nodes, a)
                return afn[t]

            # apply correction factors
            # - first and last terms have weight 1/2,
            # - even terms have opposite sign
            if len(self.nodes) % 2 == 0:
                sumA = (afn[0]/(a - self.nodes[0]) - afn[-1]/(a - self.nodes[-1])) / 2
                sumB = (1/(a - self.nodes[0]) - 1/(a - self.nodes[-1])) / 2
            else:
                sumA = (afn[0] / (a - self.nodes[0]) + afn[-1] / (a - self.nodes[-1])) / 2
                sumB = (1 / (a - self.nodes[0]) + 1 / (a - self.nodes[-1])) / 2

            # barycentric interpolation on a Chebyshev grid ...
            for tan, tfan in zip(self.nodes, afn):
                t = 1 / (a - tan)
                sumA = t * tfan - sumA
                sumB = t - sumB
            return sumA / sumB


        def derivative( self, afn ):
            '''
            return derivative of function defined by values afn
            see cheb.m in "Spectral Methods in Matlab", Ch 6, by Lloyd N. Trefethen
            use the object's barycentric differentiation matrix
            '''
            return self.diffmat @ afn


        # print string
        def __str__(self):
            '''
            user friendly print() string
            '''
            return '\'{}\' axis, min is {:.4g}, max is {:.4g}'. \
                        format(self.title, self.min_user, self.max_user )


        def __repr__(self):
            '''
            print full class data
            '''
            return 'title is \'{}\',' \
                   '\n\tuser min is {:.4g}, user max is {:.4g}' \
                   ', plot min is {:.4g}, plot max is {:.4g}' \
                   ',\n\tnodes are {}'. \
                        format(self.title, self.min_user, self.max_user, \
                               self.min, self.max, self.nodes )


    uAxis = Axis( main_params['block_params'][0]['f1_params'] )
    vAxis = Axis( main_params['block_params'][0]['f3_params'] )

    # adjust min & max for w axis
    wvals = [ func(uAxis.min_user, vAxis.min_user), func(uAxis.max_user, vAxis.min_user),
              func(uAxis.min_user, vAxis.max_user), func(uAxis.max_user, vAxis.max_user) ]
    wAxis = Axis( main_params['block_params'][0]['f2_params'], wvals )

    if "trace_init" in trace:
        print( 'lhs (u) axis: ', repr(uAxis) )
        print( 'rhs (v) axis: ', repr(vAxis) )
        print( 'middle (w) axis: ', repr(wAxis) )


    #######################################
    # call nomogram function
    # convert between plotted units and user-defined units
    # u, v and returned value are plot values,
    # call func with user-defined values
    def w(u, v):
        r = func( uAxis.user(u), vAxis.user(v) )
        return wAxis.plot(r)


    ######################################################
    #
    # the nomogram is built inside a unit square
    # pynomo scales it at the point of output.
    #
    #####################################################


    #########################################
    #
    # initialise the calculations ....
    #
    ##########################################


    w_values = np.array([[w(u, v) for v in vAxis.nodes] for u in uAxis.nodes])  # L141

    dwdu_values = np.array( [[ derivative(lambda uu: w(uu, v), u, uAxis.max, uAxis.min )
                               for v in vAxis.nodes] for u in uAxis.nodes])

    dwdv_values = np.array( [[derivative(lambda vv: w(u, vv), v, vAxis.max, vAxis.min )
                               for v in vAxis.nodes] for u in uAxis.nodes])



    #############################################################
    #
    # initial estimate of the nomogram
    #
    ##############################################################
    #
    # the u scale has uAxis.max at the top, but which way up are the w & v scales?
    #
    # draw lines from vmin to uAxis.min, uAxis.max
    # the difference in the respective w values gives the orientation of w
    # Note: if the w axis intesrsects v axis at vmin
    #       then we need to look for another point.
    #       Keep looking if that point also intersects ...
    vt = vAxis.min
    gap = 4*(vAxis.max-vAxis.min)
    while True:
        wt0 = w(uAxis.min, vt)
        wt1 = w(uAxis.max, vt)
        if not math.isclose( wt0, wt1, abs_tol=(wAxis.max-wAxis.min)/64 ):
            break
        if 'trace_init' in trace:
            print( f'{wAxis.title} axis at {wt0} touches the '
                   f'{vAxis.title} axis at {vt}, trying another point' )

        vt = vt + gap
        if vt > vAxis.max:
            if gap*32 <= vAxis.max-vAxis.min:
                sys.exit( f'too many points on {wAxis.title} '
                          f'axis near the {vAxis.title} axis - quitting' )

            gap = gap/2
            vt = vAxis.min + gap/2

    # w scale is max up (aligned to u scale) iff wt1 > wt0
    wup = wt1 > wt0
    if wup:
        yw0 = 0   # y coord of wmin
        yw1 = 1   # y coord of vAxis.max
    else:
        yw0 = 1   # y coord of wmin
        yw1 = 0   # y coord of vAxis.max


    # similar to above, draw lines from uAxis.min to vmin & vAxis.max
    # the w values for th 2 lines determine whether the v scale has
    # the same or opposite orientation to the w scale
    # look for another point if the w axis intersects u at uAxis.min, etc
    ut = uAxis.min
    gap = 4*(uAxis.max-uAxis.min)
    while True:
        wt0 = w(ut, vAxis.min)
        wt1 = w(ut, vAxis.max)
        if not math.isclose( wt0, wt1, abs_tol=(wAxis.max-wAxis.min)/64 ):
            break
        if 'trace_init' in trace:
            print( f'{wAxis.title} axis at {wt0} touches the '
                   f'{uAxis.title} axis at {ut}, trying another point' )

        ut = ut + gap
        if ut > uAxis.max:
            if gap*32 <= uAxis.max-uAxis.min:
                sys.exit( f'too many points on {wAxis.title} '
                          f'axis near the {uAxis.title} axis - quitting' )

            gap = gap/2
            ut = uAxis.min + gap/2

    vup = wup == (wt1 > wt0)

    if vup:
        # vAxis.max is at the top
        vtop = vAxis.max
        vbottom = vAxis.min
        yv0 = 0   # y coord of vmin
        yv1 = 1   # y coord of vAxis.max
    else:
        # vAxis.max is at the bottom
        vtop = vAxis.min
        vbottom = vAxis.max
        yv0 = 1   # y coord of vmin
        yv1 = 0   # y coord of vAxis.max


    ######################################################################
    # the initial condition is a nomogram that fits inside the unit square
    # - choose the u & v scales to vary linearly from (0,0) to (0,1) along
    #   the left and right edges of the unit square
    # - determine an initial estimate for the w scale

    xu = np.zeros(NN)
    yu = chebGrid
    xv = np.ones(NN)
    yv = yv0 + (yv1-yv0)*chebGrid

    # find an initial position and slope for the xw scale
    wB = w(uAxis.min,vbottom)   # bottom of w scale
    wE = w(uAxis.max,vtop)      # top of w scale
    wG = w(uAxis.max,vbottom)   # intersection of TL-BR diagonal and w wsale
    wH = w(uAxis.min,vtop)      # intersection of TR-BL diagonal and w wsale

    if "trace_init" in trace:
        print( "wup is ", wup, ", vup is ", vup,
               "wB is ", wB, ", wE is ", wE, ", wG is ", wG, ", wH is ", wH )

    if not math.isclose( wG, wH ):
        # we can use a linear estimate
        alphaxw = (wE - wG - wH + wB) / (wE - wB) / (wG - wH)
        xwB = (wH - wB) / (wE - wB) - alphaxw * (wH - wB)
        xwE = xwB + alphaxw * (wE - wB)
        if "trace_init" in trace:
            print("linear initial estimate found")
            print("xwB is ", xwB, ", xwE is ", xwE)

        # clip to unit square
        if xwE > 1:
            xwE = 1
        elif xwE < 0:
            xwE = 0
        if xwB > 1:
            xwB = 1
        elif xwB < 0:
            xwB = 0

        xw = xwB + chebGrid * (xwE - xwB)
        yw = yw0 + (yw1-yw0)*chebGrid

        if "trace_init" in trace:
            print("xw is ", xw)
            print("yw is ", yw)

    else:
        # the diagonal index lines intersect at the centre of the nomogram
        # this is a degenerate case, so use the derivatives

        # find xwE using the top isopleth, u=uAxis.max & v=vtop,
        a = (vtop-vbottom) * dwdv_values[-1][-1 if vup else 0]
        b = (uAxis.max-uAxis.min) * dwdu_values[-1][-1 if vup else 0]

        if math.isclose(a,b):
            # this also takes care of the case when a and b are both zero
            xwE = 0.5 #  put it in the middle
        else:
            xwE = a/(a+b)
        if xwE < 0:
            xwE = 0.05
        elif xwE > 1:
            xwE = 0.95

        # find wB using the bottom isopleth, u=uAxis.min & v=vbottom,
        a = (uAxis.max-uAxis.min) * dwdv_values[0][0 if vup else -1]
        b = (vtop-vbottom) * dwdu_values[0][0 if vup else -1]

        if math.isclose(a,b):
            # this also takes care of the case when a and b are both zero
            xwB = 0.5 #  put it in the middle
        else:
            xwB = a/(a+b)
        if xwB < 0:
            xwB = 0.05   # clip the lhs
        elif xwB > 1:
            xwB = 0.95   # clip the rhs

        # now fill in xw & yw
        # Notes:
        #      -  wX might not be in increasing order (ie wB might be wmax)
        #      -  use linear interpolation so xw & yw are
        #         guaranteed to lie inside the unit square

        # xw(w) goes thru (wB,xwB), (wG,0.5) & (wE,xwE)
        # yw(w) goes thru (wB,0), (wG,0.5) & (wE,1)
        if wup:
            aw = [wB, wG, wE]
            ax = [xwB, 0.5, xwE]
            ay = [0, 0.5, 1]
        else:
            aw = [wE, wG, wB]
            ax = [xwE, 0.5, xwB]
            ay = [1, 0.5, 0]
        xw = np.interp( wAxis.nodes, aw, ax )
        yw = np.interp( wAxis.nodes, aw, ay )

        if "trace_init" in trace:
            print("xw np.interp is ", xw)
            print("yw np.interp is ", yw)

        # interp1d is legacy and needs to be removed
        # for now, keep it here for checking
        if not np.allclose( xw, scipy.interpolate.interp1d( [wB, wG, wE], [xwB, 0.5, xwE],
                                         fill_value="extrapolate" )( wAxis.nodes ) ):
            sys.exit( 'xw interpolation failure' )

        if not np.allclose( yw, scipy.interpolate.interp1d( [wB, wG, wE], [0, 0.5, 1],
                                         fill_value="extrapolate" )( wAxis.nodes ) ):
            sys.exit( 'yw interpolation failure' )




    ##########################################
    #
    # this is the optimisation cost function
    #
    ##########################################

    def calc_cost(dummy):

        """
        find the cost of the candidate nomogram
        ax is a 1D array of the Chebyshev nodes of the axes

        cost is comprised of alignment accuracy, slope accuracy and area optimisation

        the accuracy component is comprised of
        - alignment error
        - gradient error

        the size component is optional
        - depends on muShape being non-zero
        - the scale of the axes
        - the angle of the index lines
        - there is a cost of exceeding the nomogram boundary

        normalise each error
        - per point pair, and
        - per range
        """

        if muShape != 0:
            lxu, lyu, lxv, lyv, lxw, lyw = np.array_split(dummy, 6)  # L657

            # the cost of straying outside the unit square:
            # the exp function gets big quickly as the points approach
            # & exceed the boundaries of the unit square.
            # apply this function to each point in the u & v axes
            # evaluate only for closest edge to each point
            # ignore cost of w axis because it is between the other 2 axes

            aa = 40 # magic scaling constant for calculating cost
                    # of approaching edge of paper area
            # the points are closest which edge?
            t1x = np.maximum( -0.3-lxu, lxu-1.3 )  # left/right edge
            t1y = np.maximum( -0.3-lyu, lyu-1.3 )  # top/bottom edge
            cost_pos_u = np.exp( aa*np.maximum(t1x, t1y) )
            t1x = np.maximum( -0.3-lxv, lxv-1.3 )
            t1y = np.maximum( -0.3-lyv, lyv-1.3 )
            cost_pos_v = np.exp( aa*np.maximum(t1x, t1y) )

        else:
            # the corner coordinates are fixed, so are not part of the optimisation
            lxu, lyu, lxv, lyv, lyw = xu, yu, xv, yv, yw
            lxu[1:-1], lyu[1:-1], lxv[1:-1], lyv[1:-1], lxw, lyw[1:-1] = \
                np.array_split( dummy, [NN-2, 2*NN-4, 3*NN-6, 4*NN-8, 5*NN-8] )  #L657

        if "trace_cost" in trace:
            print("\ntrial parameters:")
            print("lxu is ", lxu)
            print("lyu is ", lyu)
            print("lxv is ", lxv)
            print("lyv is ", lyv)
            print("lxw is ", lxw)
            print("lyw is ", lyw)
            print()

        global eAcc, eDeru, eDerv
        global eShape
        eAcc = eDeru = eDerv = eShape = 0


        # estimate position & derivative error

        fdxdu = uAxis.derivative( lxu )
        fdydu = uAxis.derivative( lyu )
        fdxdv = vAxis.derivative( lxv )
        fdydv = vAxis.derivative( lyv )
        fdxdw = wAxis.derivative( lxw )
        fdydw = wAxis.derivative( lyw )


        # precalculate coordinates and derivatives for w outside loop for speedup
        wAxis.bary.set_yi(lxw)
        txwa = wAxis.bary(w_values)

        wAxis.bary.set_yi(lyw)
        tywa = wAxis.bary(w_values)

        wAxis.bary.set_yi(fdxdw)
        dxdwa = wAxis.bary(w_values)

        wAxis.bary.set_yi(fdydw)
        dydwa = wAxis.bary(w_values)




        for iu, u in enumerate(uAxis.nodes):
            txu = lxu[iu]
            tyu = lyu[iu]
            dxdu = fdxdu[iu]
            dydu = fdydu[iu]

            for iv, v in enumerate(vAxis.nodes):
                txv = lxv[iv]
                tyv = lyv[iv]
                dxdv = fdxdv[iv]
                dydv = fdydv[iv]

                txw = txwa[iu][iv]
                tyw = tywa[iu][iv]
                dxdw = dxdwa[iu][iv]
                dydw = dydwa[iu][iv]

                tx = txu - txv
                ty = tyu - tyv
                td2 = tx * tx + ty * ty
                td = math.sqrt(td2)
                # e0 is the alignment error
                e0 = (tx * (tyu - tyw) - (txu - txw) * ty) / td  # [d]

                if td2 * resolution > 1:
                    eAcc += e0*e0                                # [d^2]

                # now determine the slope errors for the u & v axes
                tmp = dxdw * ty - tx * dydw                      # [xy]/[w]
                tuc = e0 * (tx * dxdu + ty * dydu) / td2         # [d]/[u]
                tvc = e0 * (tx * dxdv + ty * dydv) / td2         # [d]/[v]


                ### these values are predefined for speedup L141
                dwdu = dwdu_values[iu][iv]
                def quitOnNan( s, t):
                    sys.exit( """
                    error: the {} scale line around the region {} cannot be evaluated!
                    nomogen needs to evaluate the function marginally outside the given range
                    please change this limit and try again
                             """.format(s,t) )

                if math.isnan(dwdu):
                    quitOnNan(uAxis.title, uAxis.user(u))
                dedu = ((dwdu * tmp + dxdu * (tyv - tyw) - (txv - txw) * dydu)) / td - tuc  # [d]/[u] - [d]/[u]
                eDeru += dedu ** 2         # [d^2]/[u^2]


                dwdv = dwdv_values[iu][iv]
                if math.isnan(dwdv):
                    quitOnNan(vAxis.title, vAxis.user(v))
                dedv = ((dwdv * tmp + dxdv * (tyw - tyu) - (txw - txu) * dydv)) / td + tvc  # [d]/[v] + [d]/[v]
                eDerv += dedv ** 2


                check = True  # checks are slow
                check = False
                if check:
                    def gete0(au, av):
                        xau = uAxis(au, lxu)
                        yau = uAxis(au, lyu)
                        xav = vAxis(av, lxv)
                        yav = vAxis(av, lyv)
                        tw = w(au, av)
                        xtw = wAxis(tw, lxw)
                        ytw = wAxis(tw, lyw)
                        return (xau * yav + xtw * (yau - yav) + (xav - xau) * ytw - xav * yau) \
                            / math.hypot((xau - xav), (yau - yav))

                    if not math.isclose(e0, gete0(u, v), rel_tol=1e-05, abs_tol=1e-7):
                        sys.exit("gete0() error")

                    aa = v
                    tu = derivative( lambda uu: gete0(uu, aa), u, uAxis.max, uAxis.min )
                    if not math.isclose(tu, dedu, rel_tol=1e-05, abs_tol=1e-7):
                        print("(u,v) is (", u, ",", v, "), \
                        dedu is ", dedu, ", tu is ", tu, ", tuc is ", tuc)
                        sys.exit("dedu error")

                    aa = u
                    tv = derivative( lambda vv: gete0(aa, vv), v, vAxis.max, vAxis.min )
                    if not math.isclose(tv, dedv, rel_tol=1e-05, abs_tol=1e-7):
                        print("(u,v) is (", u, ",", v, "), dedv is ", dedv, ", tv is ", tv)
                        sys.exit("dedv error")


                # calculate cost of shape/measurability
                # including cost of position
                # u, v & w axes
                # use relative (%) error ie du/u
                if muShape != 0:
                    # multiply x & y coords by width & height so
                    # error is in paper units (mm), not parts of the unit square

                    # this is a factor to normalise the u variable
                    # divide by range of u to make eShape independent of u units
                    # to prevent division by zero, clip at 1% of total range
                    t = max( 0.01*(uAxis.max-uAxis.min)**2, u*u )  # [u^2]
                    xiu = (dxdu*ty - dydu*tx) * height * width     # [x y/u]
                    uShape = 1/(t*xiu**2)    # 1/([u^2] * [x^2 y^2 / u^2])

                    t = max( 0.01*(vAxis.max-vAxis.min)**2, v*v )
                    xiv = (dxdv*ty - dydv*tx) * height * width  # [x y / v]
                    vShape = 1/(t*xiv**2)   # 1/([v^2] [x^2 y^2 / v^2])

                    t = max( 0.01*(wAxis.max-wAxis.min)**2, w(u,v)**2 )
                    wShape = (t/(dwdu*xiu + dwdv*xiv))**2  # [w^2] / [w^2 x^2 y^2]

                    t = (tx * width)**2 + (ty * height)**2    # [d^2]
                    tShape = t * (uShape + vShape + wShape)   # [d^2]/[x^2 y^2]

                    tShape *= cost_pos_u[iu] + cost_pos_v[iv]
                    eShape += tShape

        # make eDeru independent of the scale of u,
        # so multiply by range of u to cancel units
        # eDeru has units (xy/u)**2
        # TODO: is there a better scaling factor?
        eDeru *= (uAxis.max - uAxis.min) ** 2                      # [d^2]
        eDerv *= (vAxis.max - vAxis.min) ** 2                      # [d^2]

        # muShape has units [d^4] so we can add
        # the shape error to the other errors.
        eShape *= muShape

        # there is one error term for each pair of points
        # so normalise to average error per point pair.
        global cur_cost
        cur_cost = (eAcc + eDeru + eDerv + eShape) / (len(uAxis.nodes) * len(vAxis.nodes))
        return cur_cost



    global itNr
    itNr = 0


    #################################################
    #
    # jacobian calculations
    #

    def calc_jac( dummy ):
        '''
        Analytically calculate the jacobian of calc_cost() above
        still a work in progress
        '''

        ## not yet complete, unused code ##

        if muShape != 0:
            lxu, lyu, lxv, lyv, lxw, lyw = np.array_split(dummy, 6)  # L657
            jxu0, jyu0, jxv0, jyv0, jxw0, jyw0 = range(0, 6*NN-1, NN)
            #print( jxu0, jyu0, jxv0, jyv0, jxw0, jyw0 )

        else:
            lxu, lyu, lxv, lyv, lyw = xu, yu, xv, yv, yw
            jxu0, jyu0, jxv0, jyv0, jxw0 = range(0, 4*NN-7, NN-2)
            jyw0 = jxw0+NN
            lxu[1:-1], lyu[1:-1], lxv[1:-1], lyv[1:-1], lxw, lyw[1:-1] = \
                np.array_split( dummy, [jyu0, jxv0, jyv0, jxw0, jyw0] )  #L657

        # share these with calc_cost()

        # precaclulate derivative functions
        fdxdu = uAxis.derivative( lxu )
        fdydu = uAxis.derivative( lyu )
        fdxdv = vAxis.derivative( lxv )
        fdydv = vAxis.derivative( lyv )
        fdxdw = wAxis.derivative( lxw )
        fdydw = wAxis.derivative( lyw )


        # precalculate coordinates and derivatives for w outside loop for speedup
        wAxis.bary.set_yi(lxw)
        txwa = wAxis.bary(w_values)
        wAxis.bary.set_yi(lyw)
        tywa = wAxis.bary(w_values)
        wAxis.bary.set_yi(fdxdw)
        dxdwa = wAxis.bary(w_values)
        wAxis.bary.set_yi(fdydw)
        dydwa = wAxis.bary(w_values)


        jac = np.zeros( len(dummy) )

        teDeru = teDerv = 0

        # accuracy error for the u scale
        jacxui = jxu0; jacyui = jyu0
        for iu in range(len(lxu)) if muShape!=0 else range(1, len(lxu)-1):
            txu = lxu[iu]
            tyu = lyu[iu]
            dxdu = fdxdu[iu]
            dydu = fdydu[iu]

            u = uAxis.nodes[iu]

            derrxui = derryui = 0   # derivative of error for xu[i]


            for iv, v in enumerate(vAxis.nodes):
                txv = lxv[iv]
                tyv = lyv[iv]
                dxdv = fdxdv[iv]
                dydv = fdydv[iv]

                txw = txwa[iu][iv]
                tyw = tywa[iu][iv]
                dxdw = dxdwa[iu][iv]
                dydw = dydwa[iu][iv]

                tx = txu - txv
                ty = tyu - tyv
                td2 = tx*tx + ty*ty
                td = math.sqrt(td2)

                tA = (txu*(tyv-tyw) + txv*(tyw-tyu) + txw*(tyu-tyv))
                tderr = 2*tA*( (tyv - tyw)*td2 - tA*tx) / (td2*td2)
                derrxui += tderr

                tderr = 2*tA*( (txw - txv)*td2 - tA*ty) / (td2*td2)
                derryui += tderr



                e0 = (tx * (tyu - tyw) - (txu - txw) * ty) / td  # [d]


                tmp = ty * dxdw - tx * dydw                      # [xy]/[w]
                tuc = e0 * (tx * dxdu + ty * dydu) / td2         # [d]/[u]
                tvc = e0 * (tx * dxdv + ty * dydv) / td2         # [d]/[v]


                ### these values are predefined for speedup L141
                dwdu = dwdu_values[iu][iv]
                quitOnNan = True               # bravely carrying on doesn't work well
                if not math.isnan(dwdu):
                    dedu = ((dwdu * tmp + (tyv - tyw) * dxdu - (txv - txw) * dydu)) / td - tuc  # [d]/[u] - [d]/[u]
                    teDeru += dedu ** 2         # [d^2]/[u^2]
                elif quitOnNan:
                    sys.exit( """
error: the {} scale line around the region u = {} cannot be represented by a finite polynomial!
       please reduce limits on this scale and try again
                                 """.format(uAxis.title, uAxis.user(u)) )

                    dwdv = dwdv_values[iu][iv]
                    if not math.isnan(dwdv):
                        dedv = ((dwdv * tmp + (tyw - tyu) * dxdv - (txw - txu) * dydv)) / td + tvc  # [d]/[v] + [d]/[v]
                        teDerv += dedv ** 2
                    elif quitOnNan:
                        sys.exit( """
error: the {} scale line around the region v = {} cannot be represented by a finite polynomial!
       try again with reduced limits on this scale
                                 """.format(vAxis.title, vAxis.user(v)) )


                    # calculate cost of shape
                    # include cost of position
                    # u & v axes
                    # divide by range of u to make eShape independent of u units, etc
                    if muShape != 0:
                        t = 0.01*(uAxis.max-uAxis.min)*(uAxis.max-uAxis.min) + u*u
                        uShape = ((dxdu*ty - dydu*tx))**2 * t   # [x^2 y^2 / u^2] * [u^2]
                        t = 0.01*(vAxis.max-vAxis.min)*(vAxis.max-vAxis.min) + v*v
                        vShape = ((dxdv*ty - dydv*tx))**2 * t
                        tShape = td2 * ( 1/uShape + 1/vShape)         # 1/[d^2]

                        #tShape *= cost_pos_u[iu] + cost_pos_v[iv]
                        #eShape += tShape

                jac[jacxui] = derrxui; jacxui += 1
                jac[jacyui] = derryui; jacyui += 1



            # accuracy error for the v scale
            jacxvi = jxv0; jacyvi = jyv0
            for iv in range(len(lxv)) if muShape!=0 else range(1, len(lxv)-1):
                txv = lxv[iv]
                tyv = lyv[iv]
                dxdv = fdxdv[iv]
                dydv = fdydv[iv]

                v = vAxis.nodes[iv]

                derrxvi = derryvi = 0   # derivative of error for xv[i]


                for iu, u in enumerate(uAxis.nodes):
                    txu = lxu[iu]
                    tyu = lyu[iu]
                    dxdu = fdxdu[iu]
                    dydu = fdydu[iu]

                    txw = txwa[iu][iv]
                    tyw = tywa[iu][iv]
                    dxdw = dxdwa[iu][iv]
                    dydw = dydwa[iu][iv]

                    tx = txu - txv
                    ty = tyu - tyv
                    td2 = tx*tx + ty*ty
                    td = math.sqrt(td2)

                    tA = (txu*(tyv-tyw) + txv*(tyw-tyu) + txw*(tyu-tyv))
                    tderr = 2*tA*( (tyw - tyu)*td2 +  tA*tx) / (td2*td2)
                    derrxvi += tderr

                    tderr = 2*tA*( (txu - txw)*td2 +  tA*ty) / (td2*td2)
                    derryvi += tderr



                    e0 = (tx * (tyu - tyw) - (txu - txw) * ty) / td  # [d]


                    tmp = ty * dxdw - tx * dydw                      # [xy]/[w]
                    tuc = e0 * (tx * dxdu + ty * dydu) / td2         # [d]/[u]
                    tvc = e0 * (tx * dxdv + ty * dydv) / td2         # [d]/[v]


                    ### these values are predefined for speedup L141
                    dwdu = dwdu_values[iu][iv]
                    quitOnNan = True               # bravely carrying on doesn't work well
                    if not math.isnan(dwdu):
                        dedu = ((dwdu * tmp + (tyv - tyw) * dxdu - (txv - txw) * dydu)) / td - tuc  # [d]/[u] - [d]/[u]
                        teDeru += dedu ** 2         # [d^2]/[u^2]
                    elif quitOnNan:
                        sys.exit( """
error: the {} scale line around the region u = {} cannot be represented by a finite polynomial!
       please reduce limits on this scale and try again
                                 """.format(uAxis.title, uAxis.user(u)) )

                    dwdv = dwdv_values[iu][iv]
                    if not math.isnan(dwdv):
                        dedv = ((dwdv * tmp + (tyw - tyu) * dxdv - (txw - txu) * dydv)) / td + tvc  # [d]/[v] + [d]/[v]
                        teDerv += dedv ** 2
                    elif quitOnNan:
                        sys.exit( """
error: the {} scale line around the region v = {} cannot be represented by a finite polynomial!
       try again with reduced limits on this scale
                                 """.format(vAxis.title, vAxis.user(v)) )


            # calculate cost of shape
            # include cost of position
            # u & v axes
            # divide by range of u to make eShape independent of u units, etc
            if muShape != 0:
                t = 0.01*(uAxis.max-uAxis.min)*(uAxis.max-uAxis.min) + u*u
                uShape = ((dxdu*ty - dydu*tx))**2 * t   # [x^2 y^2 / u^2] * [u^2]
                t = 0.01*(vAxis.max-vAxis.min)*(vAxis.max-vAxis.min) + v*v
                vShape = ((dxdv*ty - dydv*tx))**2 * t
                tShape = td2 * ( 1/uShape + 1/vShape)         # 1/[d^2]

                #tShape *= cost_pos_u[iu] + cost_pos_v[iv]
                #eShape += tShape

        jac[jacxvi] = derrxvi; jacxvi += 1
        jac[jacyvi] = derryvi; jacyvi += 1


        # now accuracy error for the w axis

        jacxwk = jxw0; jacywk = jyw0
        for kw in range(len(lxw)):  # xxx
            derrxwk = derrywk = 0   # jacobian of error for xw[k], yw[k]
            twk = wAxis.nodes[kw]
            alphak = 1
            for kk in range(len(wAxis.nodes)):  # xxx
                if kk != kw:
                    alphak = alphak*(twk - wAxis.nodes[kk])


        for iu in range(len(lxu)):  # xxx
            txu = lxu[iu]
            tyu = lyu[iu]
            dxdu = fdxdu[iu]
            dydu = fdydu[iu]

            for jv in range(len(lxv)):  # xxx

                txv = lxv[jv]
                tyv = lyv[jv]
                dxdv = fdxdv[jv]
                dydv = fdydv[jv]

                txw = txwa[iu][jv]
                tyw = tywa[iu][jv]
                dxdw = dxdwa[iu][jv]
                dydw = dydwa[iu][jv]

                tw = w_values[iu][jv]  # tw = w(u,v)

                dwk = 1/alphak
                for kk in range(len(wAxis.nodes)):
                    if kk != kw:
                        dwk = dwk*(tw - wAxis.nodes[kk])
                # dwk = d(coord at w) / d(cooord wk)

                ############## test code ###############################
                check = False
                #check = True
                if check:

                    # deltaj[j] = (-1)**j * (0.5 if j in [0, last] else 1)
                    deltaj = np.ones(len(wAxis.nodes))
                    deltaj[1::2] = -1
                    deltaj[0] = 0.5;         deltaj[-1] *= 0.5
                    #print( 'deltaj is ', deltaj )
                    if tw == twk:
                        dwk1 = 1
                    elif tw in wAxis.nodes:
                        dwk1 = 0
                    else:
                        B=0
                        for kk in range(len(wAxis.nodes)):  # xxx
                            B = B + deltaj[kk]/(tw - wAxis.nodes[kk])
                        dwk1 = (deltaj[kw]/(tw - twk))/B
                    # dwk1 = d(coord at w) / d(cooord wk)


                    A = 1; h=sqrt_eps
                    tw = tw*(1+h)
                    for kk in range(len(wAxis.nodes)):  # xxx
                        if kk != kw:
                            A = A*(tw - wAxis.nodes[kk])
                    #print( 'dwk is ', dwk, 'dwk1 is ', dwk1, ', test is ', A/alphak )

                    t = lxw[kw]
                    lxw[kw] = t+h
                    xwwt1 = wAxis(tw, lxw)
                    lxw[kw] = t-h
                    xwwt0 = wAxis(tw, lxw)
                    lxw[kw] = t
                    dddd = (xwwt1 - xwwt0) / (2*h)

                    av = (dwk + dwk1 + dddd)/3
                    if not math.isclose(av, dwk, rel_tol=1e-5, abs_tol=1e-5):
                        print( 'dwk diverges, av is ', av, 'dwk is ', dwk, ', dwk1 is ',
                               dwk1, ', dddd is ', dddd  )
                    elif not math.isclose(av, dwk1, rel_tol=1e-5, abs_tol=1e-5):
                        print( 'dwk1 diverges, av is ', av, 'dwk1 is ', dwk1  )
                    elif not math.isclose(av, dddd, rel_tol=1e-5, abs_tol=1e-5):
                        print( 'dddd diverges, av is ', av, 'dddd is ', dwk  )

                    #print( 'dwk is ', dwk, ', dddd is ', dddd, ', average is ', (dwk + dwk1 + A/B + dddd)/4 )
                ################################################

                    tx = txu - txv
                ty = tyu - tyv
                td2 = tx*tx + ty*ty

                tA = (txu*(tyv-tyw) + txv*(tyw-tyu) + txw*(tyu-tyv))
                tderr = 2*tA*dwk*(tyu - tyv) / td2
                derrxwk += tderr

                tderr = 2*tA*dwk*(txv - txu) / td2     # no y term for first & last node
                if kw not in [0, len(wAxis.nodes)-1] :
                    derrywk += tderr


        jac[jacxwk] = derrxwk; jacxwk += 1

        # first & last y nodes are fixed at 0 & 1, so they have no jacobian term
        if kw not in [0, len(wAxis.nodes)-1] :
            jac[jacywk] = derrywk; jacywk += 1



        # make eDeru independent of the scale of u,
        # so multiply by range of u to cancel units
        # eDeru has units (xy/u)**2
        teDeru *= (uAxis.max - uAxis.min) ** 2                      # [d^2]
        teDerv *= (vAxis.max - vAxis.min) ** 2                      # [d^2]

        #eShape *= muShape


        return jac / (len(uAxis.nodes) * len(vAxis.nodes))





    ##########################################


    def calc_jac00( dummy ):
        '''
        calculate the jacobian of calc_cost() above
        this version is equivalent to the current standard library function
        '''
        h = sqrt_eps                          # step size
        j = np.empty( len(dummy) )
        f0 = calc_cost( dummy )  # cost       # this has already been calculated
        for i,t in enumerate( dummy ):
            dummy[i] = t+h
            j[i] = ( calc_cost(dummy)-f0 ) / h
            dummy[i] = t

        return j


    def calc_jac2( dummy ):
        '''
        calculate the jacobian of calc_cost() above
        more accurate version (O(h^2)) than calc_jac00(), but
        more expensive, requires 2 evaluations of the cost function per loop
        seems equivalent to jac='3-point'
        '''
        h = sqrt_eps*4                             # step size
        j = np.empty( len(dummy) )
        for i,t in enumerate( dummy ):
            dummy[i] = t-h
            f0 = calc_cost( dummy )
            dummy[i] = t+h
            j[i] = ( calc_cost(dummy)-f0 ) / (2*h)
            dummy[i] = t

        if False and itNr % 10 == 0:
            print( '\njac2 is ', j )
            print( 'jac00 is ', calc_jac00( dummy ) )
            print( 'jac1 is ', calc_jac( dummy ), '\n' )

        return j


    ##################################################################
    #
    # run this function whenever the optimise function enters a
    # new iteration

    def newStep(xk):
        global itNr
        global cur_cost, old_cost
        global lcost1

        itNr += 1
        print("\r", main_params['filename'], ": {:3d} ".format(itNr), sep="", end="")

        # display the details if there is any tracing
        if trace:
            print(", eAcc {:.2e}".format(eAcc), end='')
            print(", eDer {:.2e}, {:.2e}".format(eDeru, eDerv), end='')
            if muShape != 0:
                print(", eShape {:.2e}".format(eShape), end='')
            if old_cost != 0 and cur_cost < old_cost:
                print(", cost improvement is {:2.0f}%".format(100 * (old_cost - cur_cost) / old_cost), end='')
            old_cost = cur_cost
        else:
            # draw a progress bar
            # assume final cost is 5e-7, => ln(5e-7) is  -14.5
            # draw bar in log proportion to final cost
            wid = 70  # width of bar
            if itNr == 1:
                p = 0   # start proportion
                lcost1 = math.log(cur_cost)  # start cost
            else:
                # calculate the proportion for current cost
                p = round( wid*(math.log(cur_cost) - lcost1) / (-14.5 - lcost1) )
            print( '#'*p, '-'*(wid-p), '|', sep='', end='' )


    #
    # the nomogram is given by the coordinates of the node points that
    # minimises the error
    #

    assert len(xu) == NN
    if muShape != 0:
        x0 = np.concatenate([xu, yu, xv, yv, xw, yw])  # L657
    else:
        x0 = np.concatenate([xu[1:-1], yu[1:-1],
                             xv[1:-1], yv[1:-1],
                             xw, yw[1:-1]])  # L657


    # BFGS is the default optimisation method - it is most reliable, but might be slow
    # TODO: consider non linear least squares method to speed up
    # increase gtol if this terminates due to loss of precision
    res = scipy.optimize.minimize( calc_cost, x0,
                                   method='BFGS', callback=newStep, jac=calc_jac2,
                                   options={'disp': False, 'gtol': 1e-4, 'maxiter': None},
                                   tol=1e-4 )


    ######################################################
    #
    # the nomogram is solved - now check & report results
    #
    #######################################################

    print()
    if "trace_result" in trace:
        print(res)
    else:
        print( res.nit, " iterations, cost function is {:.2e}, ".format(res.fun), res.message, sep='')

    if res.success or 'precision loss' in res.message:  # lack of precision still gives a result
        if muShape != 0:
            xu, yu, xv, yv, xw, yw = \
                np.array_split(res.x, 6)  # L657
        else:
            xu[1:-1], yu[1:-1], xv[1:-1], yv[1:-1], xw, yw[1:-1] = \
                np.array_split(res.x,
                               [NN - 2, 2 * NN - 4, 3 * NN - 6,
                                4 * NN - 8, 5 * NN - 8])  # L657
        #print( '\njacobian is:\n', res.jac, '\n' )
        #print( '\ncalc_jac0 is:\n', calc_jac0(res.x), '\n' )

        if "trace_result" in trace:
            print("solution is ...")
            print("xu is ", xu)
            print("yu is ", yu)
            print("xv is ", xv)
            print("yv is ", yv)
            print("xw is ", xw)
            print("yw is ", yw)


        # check derivatives at each node pair...
        fdxdu = uAxis.derivative( xu )
        fdydu = uAxis.derivative( yu )
        fdxdv = vAxis.derivative( xv )
        fdydv = vAxis.derivative( yv )
        fdxdw = wAxis.derivative( xw )
        fdydw = wAxis.derivative( yw )

        maxerr = 0

        wAxis.bary.set_yi(fdxdw)
        dxdwa = wAxis.bary(w_values)
        wAxis.bary.set_yi(fdydw)
        dydwa = wAxis.bary(w_values)
        wAxis.bary.set_yi(xw)
        xwcoorda = wAxis.bary(w_values)
        wAxis.bary.set_yi(yw)
        ywcoorda = wAxis.bary(w_values)

        for i, (xucoord, yucoord, dxdu, dydu) in enumerate( zip(xu, yu, fdxdu, fdydu) ):

            for j, (xvcoord, yvcoord, dxdv,  dydv) in enumerate( zip(xv, yv, fdxdv,fdydv) ):
                wvalue = w_values[i,j]
                xwcoord = xwcoorda[i,j]
                ywcoord = ywcoorda[i,j]
                dxdw = dxdwa[i,j]
                dydw = dydwa[i,j]

                dwdu = dwdu_values[i][j]
                dwdv = dwdv_values[i][j]
                lhs = dwdu * ((xucoord - xvcoord) * dydw - (yucoord - yvcoord) * dxdw)
                rhs = (xwcoord - xvcoord) * dydu - (ywcoord - yvcoord) * dxdu
                t = abs(lhs - rhs)*(uAxis.max-uAxis.min)
                maxerr = max(maxerr, t)

                lhs = dwdv * ((xucoord - xvcoord) * dydw - (yucoord - yvcoord) * dxdw)
                rhs = (xucoord - xwcoord) * dydv - (yucoord - ywcoord) * dxdv
                t = abs(lhs - rhs)*(vAxis.max-vAxis.min)
                maxerr = max(maxerr, t)

                lhs = dwdv * ((xwcoord - xvcoord) * dydu - (ywcoord - yvcoord) * dxdu)
                rhs = dwdu * ((xucoord - xwcoord) * dydv - (yucoord - ywcoord) * dxdv)
                t = abs(lhs - rhs)*(uAxis.max-uAxis.min)*(vAxis.max-vAxis.min)/(wAxis.max-wAxis.min)
                maxerr = max(maxerr, t)
        print("max derivative error is {:.2g}".format(maxerr))


        # now iterate over the length of the u & v scales,
        # verifying the nomogram at each step

        hmm = math.sqrt(height * width) # for converting distance from fractions of unit square to mm

        # verify alignment every d mm
        d = 1  # choose 1mm
        print("checking solution every ", d, "mm, ")
        ds = d / hmm  # step size as fraction of the unit square

        upoints = []
        u = uAxis.min
        while u <= uAxis.max:
            #print("u is ", u, 100*(u-uAxis.min)/(uAxis.max-uAxis.min), "%")
            upoints.append(u)
            dx = uAxis( u, fdxdu )
            dy = uAxis( u, fdydu )
            # next u:
            dd = ds / math.hypot(dx, dy)
            if uAxis.max-dd < u < uAxis.max:
                u = uAxis.max
            else:
                u = u+dd

        vpoints = []
        v = vAxis.min
        while v <= vAxis.max:
            #print("v is ", v, 100*(v-vAxis.min)/(vAxis.max-vAxis.min), "%")
            vpoints.append(v)
            dx = vAxis(v, fdxdv )
            dy = vAxis(v, fdydv )
            # next v:
            dd = ds / math.hypot(dx, dy)
            if vAxis.max - dd < v < vAxis.max:
                v = vAxis.max
            else:
                v = v+dd


        uAxis.bary.set_yi(xu)
        xucoorda = uAxis.bary(upoints)
        uAxis.bary.set_yi(yu)
        yucoorda = uAxis.bary(upoints)

        vAxis.bary.set_yi(xv)
        xvcoorda = vAxis.bary(vpoints)
        vAxis.bary.set_yi(yv)
        yvcoorda = vAxis.bary(vpoints)

        logAlignmentErrors = 'LogAlignment' in main_params['block_params'][0] \
            and main_params['block_params'][0]['LogAlignment']
        if logAlignmentErrors:
            fn = main_params['filename']+'.error.py'
            print( 'writing alignment errors to {}'.format(fn) )
            try:
                logFile = open( fn, 'w')
            except (PermissionError, OSError) as e:
                print("Cannot open error log file: '{}'".format(e))
                print('Alignment log not saved')
                logAlignmentErrors = False
            else:
                logFile.write( '# nomogen alignment error key value pairs, etc.  error is in mm\n' )
                logFile.write( '# generated by nomogen {} \n'.
                               format( datetime.datetime.now().strftime('%a %d %b %Y %H:%M:%S'), ))
                logFile.write( "nomo_id = {{ 'file':'{}', 'title': '{}', 'paper_height': {}, 'paper_width':{} }}\n". \
                               format( main_params['filename'], wAxis.title,  # ???
                                       main_params['paper_height'],
                                       main_params['paper_width'], ) )
                logFile.write( 'Nomo_error = [\n')

        maxdiff = 0
        for u, xucoord, yucoord in zip( upoints, xucoorda, yucoorda ):
            for v, xvcoord, yvcoord in zip( vpoints, xvcoorda, yvcoorda ):
                def report_alignment(difference, xwcoord, ywcoord):
                    print("u is ", uAxis.user(u), " at pos (", xucoord, ",", yucoord, ")")
                    print("v is ", vAxis.user(v), " at pos (", xvcoord, ",", yvcoord, ")")
                    print("w is ", wAxis.user(wvalue), " at pos (", xwcoord, ",", ywcoord, ")")
                    print( "alignment difference is {:5.2g} about {:5.2f} mm".format(difference, difference * hmm) )

                wvalue = w(u, v)
                if wAxis.params['u_min'] < wAxis.user(wvalue) < wAxis.params['u_max']:
                    xwcoord = wAxis(wvalue, xw)
                    ywcoord = wAxis(wvalue, yw)
                    difference = ((xucoord - xvcoord) * (yucoord - ywcoord) -
                                  (xucoord - xwcoord) * (yucoord - yvcoord)) / \
                                  math.hypot(xucoord - xvcoord, xvcoord - yvcoord)

                    if logAlignmentErrors:
                        logFile.write( "{{'w': {}, 'error': {}, 'u': {}, 'v': {}}},\n".
                                       format(wAxis.user(wvalue), hmm*difference, u, v) )

                    difference = abs(difference)
                    if wvalue < wAxis.min and not math.isclose(wvalue, wAxis.min, rel_tol=0.01):
                        report_alignment(difference, xwcoord, ywcoord)
                        print( "scale range error, please check w scale min limits",
                               "\nw({:g}, {:g}) = {:g} < wAxis.min, {:g}".
                               format(uAxis.user(u), vAxis.user(v), wAxis.user(wvalue), wAxis.min_user))
                        #sys.exit("scale range error, please check w scale min limits"
                            #     "\nw({:g}, {:g}) = {} < wAxis.min, {}".format(u, v, wvalue, wAxis.min))
                    elif wvalue > wAxis.max and not math.isclose(wvalue, wAxis.max, rel_tol=0.01):
                        report_alignment(difference, xwcoord, ywcoord)
                        print("scale range error, please check w scale max limits",
                              "\nw({:g}, {:g}) = {:g} > wAxis.max, {:g}".
                              format(uAxis.user(u), vAxis.user(v), wAxis.user(wvalue), wAxis.max_user))
                        #sys.exit("scale range error, please check w scale max limits"
                        #         "\nw({:g}, {:g}) = {} > wAxis.max, {}".format(u, v, wvalue, wAxis.max))

                    if difference > maxdiff:
                        maxdiff = difference
                        if "trace_alignment" in trace:
                            report_alignment(difference, xwcoord, ywcoord)

        if logAlignmentErrors:
            logFile.write( ']\n')
            logFile.close()


        aler = maxdiff * math.sqrt(width * height)
        print("alignment error is estimated at less than {:5.2g} mm".format(aler))
        if aler > 0.2:
            print("alignment errors are possible - please check.")
            print( "This nomogram used a polynomial defined with {} points ".format(NN) )
            print("Try increasing this, or reduce the range of one or more scales")

        # return the resulting axes into the nomogram parameters
        uAxis.params.update({'f': lambda u: uAxis(uAxis.plot(u), xu),
                             'g': lambda u: uAxis(uAxis.plot(u), yu),
                             'h': lambda u: 1.0,
                             })
        vAxis.params.update({'f': lambda v: vAxis(vAxis.plot(v), xv),
                             'g': lambda v: vAxis(vAxis.plot(v), yv),
                             'h': lambda v: 1.0,
                             })
        wAxis.params.update({'f': lambda w: wAxis(wAxis.plot(w), xw),
                             'g': lambda w: wAxis(wAxis.plot(w), yw),
                             'h': lambda w: 1.0,
                             })

        # get & print footer info
        if 'footer_string' in main_params:
            txt = r'\tiny \hfil {} \hfil'.format( main_params['footer_string'])
        else:
            datestr = datetime.datetime.now().strftime("%d %b %y")
            if aler > 0.1:
                tolstr = r", est tolerance {:5.2g} mm".format(aler)
            else:
                tolstr = ""

            # '\' char is escape, '_' char is subscript, '$' is math mode, etc
            escapes = "".maketrans({ '\\': r'$ \backslash $',
                                     '^': r'\^{}',
                                     '_': r'\_',
                                     '~': r'$ \sim $',
                                     '$': r'\$',
                                     '#': r'\#',
                                     '%': r'\%',
                                     '&': r'\&',
                                     '{': r'\{',
                                     '}': r'\}' })
            txt = r'\tiny \hfil {}: created by nomogen {} {} \hfil'. \
                  format( main_params['filename'].translate(escapes), datestr, tolstr)
        #print("txt is \"", txt, "\"", sep='')

        footerText = {'x': 0,
                      'y': 0.0,
                      'nomogen_footer' : True,
                      'text': txt,
                      'width': width/10,
                     }
        if 'extra_texts' not in main_params:
            main_params['extra_texts'] = [footerText]
        else:
            # look for an existing nomogen footer in the extra texts
            # replace if found, otherwise add
            for i,xt in enumerate( main_params['extra_texts'] ):
                if 'nomogen_footer' in xt:
                    main_params['extra_texts'][i] = footerText
                    break
            else:
                main_params['extra_texts'].append(footerText)


        ##################################################
        #
        # check if the tic marks on the scales clash
        # swap sides if the w scale is too close to one of the outer scales
        # unless the user has already specified a preference ...

        # account for scales increasing downwards
        #print( f'wAxis.min_user is {wAxis.min_user}, wAxis.max_user is {wAxis.max_user}, wAxis.plot(wAxis.max_user)is {wAxis.plot(wAxis.max_user)}' )
        if (yw[0] > yw[-1]) != (wAxis.plot(wAxis.min_user) > wAxis.plot(wAxis.max_user)):
            txw = xw[::-1]
            wEast = 'left'
            wWest = 'right'
        else:
            txw = xw
            wEast = 'right'
            wWest = 'left'

        if (yv[0] > yv[-1]) != (vAxis.plot(vAxis.min_user) > vAxis.plot(vAxis.max_user)):
            txv = xv[::-1]
            vEast = 'left'
        else:
            txv = xv
            vEast = 'right'

        if uAxis.plot(uAxis.min_user) > uAxis.plot(uAxis.max_user):
            #print( 'u scale increases downward' )
            txu = xu[::-1]
            uWest = 'right'
        else:
            #print( 'u scale increases upward' )
            txu = xu
            uWest = 'left'

        # the distance between axes
        distuw = txw - txu
        distvw = txv - txw
        #print( 'distance u..w is {}, distance v..w is {}'.format(distuw, distvw) )

        if ('tick_side' not in uAxis.params) and (np.min(distuw) < 0.2):
            print( 'putting u scale ticks on {} side'.format(uWest) )
            uAxis.params.update({'tick_side': uWest})
            uAxis.params.update({'turn_relative': True})

        if ('tick_side' not in vAxis.params) and (np.min(distvw) < 0.2):
            print( 'putting v scale ticks on {} side'.format(vEast) )
            vAxis.params.update({'tick_side': vEast})
            vAxis.params.update({'turn_relative': True})


        #print( 'vw norm is ', np.linalg.norm(xv - xw), ', uw norm is ', np.linalg.norm(xw - xu) )
        if not 'tick_side' in wAxis.params:
            if np.min(distuw) < 0.2:
                wAxis.params.update({'tick_side': wEast})
            elif np.min(distvw) < 0.2:
                wAxis.params.update({'tick_side': wWest})
            elif np.linalg.norm(distuw) > np.linalg.norm(distvw):
                # the w scale line is closer to the v scale line
                wAxis.params.update({'tick_side': wWest})
            else:
                # the w scale line is closer to the u scale line
                wAxis.params.update({'tick_side': wEast})
            wAxis.params.update({'turn_relative': True})
            print( 'putting w scale ticks on {} side'.format(wAxis.params['tick_side']) )


        ########################################################
        #
        # check for dual scales
        # the first block holds the type_9 axes, (uAxis.params, etc)
        # then type_8 or maybe type_9 blocks follow with the dual scales

        ## TODO:
        # atm the dual scales are considered after calculating the axes.
        # in a future version, more than a single equation might be allowed,
        # so the blocks should be analysed before the scales are calculated

        for params in [uAxis.params, vAxis.params, wAxis.params]:
            if 'tag' in params:
                # look for a matching tag
                # check and copy the newly calculated functions across
                ltag = params['tag']
                for b in main_params['block_params'][1:]:
                    if b['block_type'] == 'type_8':
                        if not 'f_params' in b:
                            sys.exit('''
                            dual axis tag '{}' should have an 'f_params' parameter
                            '''.format(ltag) )
                        laxis = b['f_params']
                        if 'tag' in laxis and laxis['tag'] == ltag:
                            # don't copy function if the user has already defined one
                            if not 'function_x' in laxis:
                                print( 'dual scales found for axes with tag \'{}\''.format(ltag) )
                                if not 'align_func' in laxis:
                                    sys.exit('''
                                    dual axis tag '{}' must have an 'align_func' parameter
                                    '''.format(ltag) )
                                fal = laxis['align_func']
                                # params & fal must be evaluated now, not when 'fuction_?' is called
                                laxis['function_x'] = functools.partial(lambda u, p, f: p['f'](f(u)),
                                                                        p=params, f=fal)
                                laxis['function_y'] = functools.partial(lambda u, p, f: p['g'](f(u)),
                                                                        p=params, f=fal)

                    elif b['block_type'] == 'type_9':
                        for fp in [ 'f1_params', 'f2_params', 'f3_params' ]:
                            if not fp in b:
                                sys.exit('''
                                dual axis tag '{}' should have an '{}' parameter
                                '''.format(ltag, fp) )
                            laxis = b[fp]
                            # don't copy function if the user has already defined one
                            if not 'f' in laxis:
                                if 'tag' in laxis and laxis['tag'] == ltag:
                                    print( 'dual scales found for axes with tag \'{}\''.format(ltag) )
                                    if not 'align_func' in laxis:
                                        sys.exit('''
                                        dual axis tag '{}' must have an 'align_func' parameter
                                        '''.format(ltag) )
                                    fal = laxis['align_func']
                                    # params & fal must be evaluated now, not when 'fuction_?' is called
                                    laxis['f'] = functools.partial(lambda u, p, f: p['f'](f(u)),
                                                                            p=params, f=fal)
                                    laxis['g'] = functools.partial(lambda u, p, f: p['g'](f(u)),
                                                                            p=params, f=fal)
                                    laxis['h'] = params['h']


# end of Nomogen()


######################### end of nomogen.py ########################
