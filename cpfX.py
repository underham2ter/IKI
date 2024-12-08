"""
Rational Approximations for the Complex Error Function (aka Complex Probability Function cpf)
w(z)  =  (i/pi) integral exp(-t**2) / (z-t) dt  =  K(x,y) + iL(x,y)

Josef Humlicek:
An efficient method for evaluation of the complex probability function: the Voigt function and its derivatives.
J. Quant. Spectrosc. & Radiat. Transfer, 21, 309-313, 1979;  doi: 10.1016/0022-4073(79)90062-1

J. Humlicek:
Optimized computation of the Voigt and complex probability function.
J. Quant. Spectrosc. & Radiat. Transfer, 27, 437­444, 1982;  doi: 10.1016/0022-4073(82)90078-4

This file is part of the supplementary of
F. Schreier:
The Voigt and complex error function: Humlicek's rational approximation generalized.
MNRAS 479, 3068–3075 (2018), doi: 10.1093/mnras/sty1680

This file includes:
 * Python/NumPy implementations of the original Fortran 77 source code
   (with two versions to combine the two regions: where and vectorize)
 * a generalization of the region I rational approximation to arbitrary (even) number of terms
 * an optimized implementation of the region I approximation for 16 terms using a single fraction
   The a and b coefficients of the numerator, denominator polynomials have been evaluated using sympy

Most functions here use a single complex argument z=x+i*y    (except for the last two: hum1zpf16m, hum2zpf16m)
"""


####################################################################################################################################
#####                                                                                                                          #####
#####    Changes November 2022                                                                                                 #####
#####    * comment header replaced by docstring                                                                                #####
#####    * __all__ statement to prevent import of constants                                                                    #####
#####    * loop version of zpf16                                                                                               #####
#####    * timeit benchmarks for various versions of the zpf16 approximation (bottom of file)                                  #####
#####                                                                                                                          #####
####################################################################################################################################

import numpy as np

try:  from scipy.special import roots_hermite  # added recently, old versions do not have it
except ImportError:  print ('WARNING:  "from scipy.special import roots_hermite" failed, no cpfX function!')

# public functions (i.e. do not import the constants):
__all__ = 'cpf12v cpf12w  cpfX  zpf16p zpf16h zpf16l zpf16P  hum1zpf16w hum2zpf16w hum1zpf16m hum2zpf16m'.split()


recPi = 1/np.pi
recSqrtPi  = np.sqrt(recPi)
iRecSqrtPi = 1j*recSqrtPi

###################################################################################################


#################################

cr=1.5; crr=2.25  #  y0 and y0**2 of the original (y0 corresponds to delta below)
ct = np.array([0., .3142403762544,  .9477883912402, 1.5976826351526, 2.2795070805011, 3.0206370251209, 3.88972489786978])
ca = np.array([0., -1.393236997981977,    -0.2311524061886763,   +0.1553514656420944,
                   -0.006218366236965554, -9.190829861057117e-5, +6.275259577e-7])
cb = np.array([0.,  1.011728045548831,    -0.7519714696746353, 0.01255772699323164,
                    0.01002200814515897, -2.420681348155727e-4,  5.008480613664576e-7])


def cpf12a (z):
	""" Humlicek (1979) complex probability function:  rational approximation for y>0.85 OR |x|<18.1*y+1.65  (region I).
	    Eq. (6)
	"""
	x, y = z.real, z.imag
	ry = cr+y
	ryry = ry**2
	wk =  (ca[1]*(x-ct[1]) + cb[1]*ry) / ((x-ct[1])**2 + ryry) - (ca[1]*(x+ct[1]) - cb[1]*ry) / ((x+ct[1])**2 + ryry) \
	    + (ca[2]*(x-ct[2]) + cb[2]*ry) / ((x-ct[2])**2 + ryry) - (ca[2]*(x+ct[2]) - cb[2]*ry) / ((x+ct[2])**2 + ryry) \
	    + (ca[3]*(x-ct[3]) + cb[3]*ry) / ((x-ct[3])**2 + ryry) - (ca[3]*(x+ct[3]) - cb[3]*ry) / ((x+ct[3])**2 + ryry) \
	    + (ca[4]*(x-ct[4]) + cb[4]*ry) / ((x-ct[4])**2 + ryry) - (ca[4]*(x+ct[4]) - cb[4]*ry) / ((x+ct[4])**2 + ryry) \
	    + (ca[5]*(x-ct[5]) + cb[5]*ry) / ((x-ct[5])**2 + ryry) - (ca[5]*(x+ct[5]) - cb[5]*ry) / ((x+ct[5])**2 + ryry) \
	    + (ca[6]*(x-ct[6]) + cb[6]*ry) / ((x-ct[6])**2 + ryry) - (ca[6]*(x+ct[6]) - cb[6]*ry) / ((x+ct[6])**2 + ryry)
	wl =  (cb[1]*(x-ct[1]) - ca[1]*ry) / ((x-ct[1])**2 + ryry) + (cb[1]*(x+ct[1]) + ca[1]*ry) / ((x+ct[1])**2 + ryry) \
	    + (cb[2]*(x-ct[2]) - ca[2]*ry) / ((x-ct[2])**2 + ryry) + (cb[2]*(x+ct[2]) + ca[2]*ry) / ((x+ct[2])**2 + ryry) \
	    + (cb[3]*(x-ct[3]) - ca[3]*ry) / ((x-ct[3])**2 + ryry) + (cb[3]*(x+ct[3]) + ca[3]*ry) / ((x+ct[3])**2 + ryry) \
	    + (cb[4]*(x-ct[4]) - ca[4]*ry) / ((x-ct[4])**2 + ryry) + (cb[4]*(x+ct[4]) + ca[4]*ry) / ((x+ct[4])**2 + ryry) \
	    + (cb[5]*(x-ct[5]) - ca[5]*ry) / ((x-ct[5])**2 + ryry) + (cb[5]*(x+ct[5]) + ca[5]*ry) / ((x+ct[5])**2 + ryry) \
	    + (cb[6]*(x-ct[6]) - ca[6]*ry) / ((x-ct[6])**2 + ryry) + (cb[6]*(x+ct[6]) + ca[6]*ry) / ((x+ct[6])**2 + ryry)
	return wk+1j*wl   # wk, wl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def cpf12b (z):
	""" Humlicek (1979) complex probability function:  rational approximation for y<0.85 AND |x|>18.1*y+1.65  (region II).
	    Eq. (11)
	"""
	x, y = z.real, z.imag
	ry   = cr+y      # yy0   = y+1.5
	y2r  = y +2.*cr  # y2y0  = y+3.0
	rry  = cr*ry     # y0yy0 = 1.5*(y+1.5)
	ryry = ry**2     # yy0sq = (y+1.5)**2
	wk =  ( cb[1]*((x-ct[1])**2-rry) - ca[1]*(x-ct[1])*y2r ) / (((x-ct[1])**2+crr)*((x-ct[1])**2+ryry)) \
            + ( cb[1]*((x+ct[1])**2-rry) + ca[1]*(x+ct[1])*y2r ) / (((x+ct[1])**2+crr)*((x+ct[1])**2+ryry)) \
	    + ( cb[2]*((x-ct[2])**2-rry) - ca[2]*(x-ct[2])*y2r ) / (((x-ct[2])**2+crr)*((x-ct[2])**2+ryry)) \
            + ( cb[2]*((x+ct[2])**2-rry) + ca[2]*(x+ct[2])*y2r ) / (((x+ct[2])**2+crr)*((x+ct[2])**2+ryry)) \
	    + ( cb[3]*((x-ct[3])**2-rry) - ca[3]*(x-ct[3])*y2r ) / (((x-ct[3])**2+crr)*((x-ct[3])**2+ryry)) \
            + ( cb[3]*((x+ct[3])**2-rry) + ca[3]*(x+ct[3])*y2r ) / (((x+ct[3])**2+crr)*((x+ct[3])**2+ryry)) \
	    + ( cb[4]*((x-ct[4])**2-rry) - ca[4]*(x-ct[4])*y2r ) / (((x-ct[4])**2+crr)*((x-ct[4])**2+ryry)) \
            + ( cb[4]*((x+ct[4])**2-rry) + ca[4]*(x+ct[4])*y2r ) / (((x+ct[4])**2+crr)*((x+ct[4])**2+ryry)) \
	    + ( cb[5]*((x-ct[5])**2-rry) - ca[5]*(x-ct[5])*y2r ) / (((x-ct[5])**2+crr)*((x-ct[5])**2+ryry)) \
            + ( cb[5]*((x+ct[5])**2-rry) + ca[5]*(x+ct[5])*y2r ) / (((x+ct[5])**2+crr)*((x+ct[5])**2+ryry)) \
	    + ( cb[6]*((x-ct[6])**2-rry) - ca[6]*(x-ct[6])*y2r ) / (((x-ct[6])**2+crr)*((x-ct[6])**2+ryry)) \
            + ( cb[6]*((x+ct[6])**2-rry) + ca[6]*(x+ct[6])*y2r ) / (((x+ct[6])**2+crr)*((x+ct[6])**2+ryry))
	wl =  (cb[1]*(x-ct[1]) - ca[1]*ry) / ((x-ct[1])**2 + ryry) + (cb[1]*(x+ct[1]) + ca[1]*ry) / ((x+ct[1])**2 + ryry) \
	    + (cb[2]*(x-ct[2]) - ca[2]*ry) / ((x-ct[2])**2 + ryry) + (cb[2]*(x+ct[2]) + ca[2]*ry) / ((x+ct[2])**2 + ryry) \
	    + (cb[3]*(x-ct[3]) - ca[3]*ry) / ((x-ct[3])**2 + ryry) + (cb[3]*(x+ct[3]) + ca[3]*ry) / ((x+ct[3])**2 + ryry) \
	    + (cb[4]*(x-ct[4]) - ca[4]*ry) / ((x-ct[4])**2 + ryry) + (cb[4]*(x+ct[4]) + ca[4]*ry) / ((x+ct[4])**2 + ryry) \
	    + (cb[5]*(x-ct[5]) - ca[5]*ry) / ((x-ct[5])**2 + ryry) + (cb[5]*(x+ct[5]) + ca[5]*ry) / ((x+ct[5])**2 + ryry) \
	    + (cb[6]*(x-ct[6]) - ca[6]*ry) / ((x-ct[6])**2 + ryry) + (cb[6]*(x+ct[6]) + ca[6]*ry) / ((x+ct[6])**2 + ryry)
	return np.exp(-x*x)+y*wk+1j*wl   # wk, wl


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def _cpf12v_ (z):
	"""  Combination of Humlicek's (1979) rational approximations (scalar argument z only, for internal use). """
	if z.imag>0.85 or abs(z.real)<18.1*z.imag+1.65:  return cpf12a(z)
	else:                                            return cpf12b(z)  # speed roughly identical for np.exp and math.exp


cpf12v = np.vectorize(_cpf12v_)                                            # factor 400 slower than cpf12w


def cpf12w (z):
	""" Humlicek (JQSRT 1979): An efficient method for evaluation of the complex probability function
	    Combination of the region I and II rational approximations using numpy.where
	"""
	xMask = abs(z.real)<18.1*z.imag+1.65
	yMask = z.imag>0.85
	mask  = xMask*yMask
	return np.where (mask, cpf12a(z), cpf12b(z))


####################################################################################################################################
#####    Generalization of Humlicek's region I rational approximation                                                          #####
####################################################################################################################################

def cpfX (z, n=12, delta=1.5):
	""" Humlicek (1979) complex probability function:  a rational approximation for all z=x+iy
	    Eq. (6)
	    WARNING:  large errors for small y and n<16
	"""
	if n%2>0:  raise SystemExit ("ERROR --- cpfX:  n is odd")
	x, y = z.real, z.imag
	roots, weights = roots_hermite(n)
	alfa = -recPi * weights * np.exp(delta**2) * np.sin(2*roots*delta)
	beta =  recPi * weights * np.exp(delta**2) * np.cos(2*roots*delta)
	ry = delta+y
	ryry = ry**2

	wk=0.0;  wl=0.0
	for t,a,b in zip(roots[n//2:],alfa[n//2:],beta[n//2:]):         # old Python2:  zip(roots[n/2:],alfa[n/2:],beta[n/2:]):
		wk += (a*(x-t) + b*ry) / ((x-t)**2 + ryry) - (a*(x+t) - b*ry) / ((x+t)**2 + ryry)
		wl += (b*(x-t) - a*ry) / ((x-t)**2 + ryry) + (b*(x+t) + a*ry) / ((x+t)**2 + ryry)
	return wk+1j*wl   # wk, wl


####################################################################################################################################
#####       Optimized (single fraction) Humlicek region I rational approximation for n=16 and delta=1.31183                    #####
####################################################################################################################################

aa  = np.array([ 41445.0374210222,
                 -136631.072925829j,
                 -191726.143960199,
                 268628.568621291j,
                 173247.907201704,
                 -179862.56759178j,
                 -63310.0020563537,
                 56893.7798630723j,
                 11256.4939105413,
                 -9362.62673144278j,
                 -1018.67334277366,
                 810.629101627698j,
                 44.5707404545965,
                 -34.5401929182016j,
                 -0.740120821385939,
                 0.564189583547714j])  # identical to 1/sqrt(pi) except for last two digits

bb  = np.array([ 7918.06640624997, 0.0,
                 -126689.0625,     0.0,
                 295607.8125,      0.0,
                 -236486.25,       0.0,
                 84459.375,        0.0,
                 -15015.0,         0.0,
                 1365.0,           0.0,
                 -60.0,            0.0,
                 1.0])


# The poly1d class assumes that the very first coefficient corresponds to the highest power of z and the last coefficient to z**0
numPoly16 = np.poly1d(np.flipud(aa))
denPoly16 = np.poly1d(np.flipud(bb))        # %timeit denPoly16(z)  10000 loops, best of 3:  28 microseconds per loop
denPoly8e = np.poly1d(np.flipud(bb[::2]))   # %timeit denPoly8(z*z) 100000 loops, best of 3: 16.5 microseconds per loop

##### for time benchmarks see bottom of file

def zpf16p (z):
	""" Humlicek (1979) complex probability function for n=16 and delta=1.31183.
	    Optimized rational approximation using numpy.poly1d  (applicable for all z).

	    Maximum relative error to wofz in 0<x<25 and 1e-8<y<1e2:  7.747e-5 and 4.08e-6 for K and L
	    """
	Z = z + 1.31183j
	ZZ = Z*Z
	return numPoly16(Z)/denPoly8e(ZZ)    # %timeit zpf16p    10000 loops, best of 3: 173 microseconds per loop


def zpf16h (z):
	""" Humlicek (1979) complex probability function for n=16 and delta=1.31183.
	    Optimized rational approximation with Horner scheme  (applicable for all z). """
	Z = z + 1.31183j
	ZZ = Z*Z
	return (((((((((((((((aa[15]*Z+aa[14])*Z+aa[13])*Z+aa[12])*Z+aa[11])*Z+aa[10])*Z+aa[9])*Z+aa[8])*Z+aa[7])*Z+aa[6])*Z+aa[5])*Z+aa[4])*Z+aa[3])*Z+aa[2])*Z+aa[1])*Z+aa[0]) \
	      / ((((((((ZZ+bb[14])*ZZ+bb[12])*ZZ+bb[10])*ZZ+bb[8])*ZZ+bb[6])*ZZ+bb[4])*ZZ+bb[2])*ZZ+bb[0])


def zpf16l (z):
	""" Humlicek (1979) complex probability function for n=16 and delta=1.31183.
	    Optimized rational approximation with Horner scheme rewritten as a for loop (applicable for all z). """
	Z = z + 1.31183j
	ZZ = Z*Z
	numer, denom = aa[-1], 0.0
	for a in aa[-2::-1]:  numer = numer * Z + a
	for b in bb[::-2]:    denom = denom * ZZ + b
	return numer / denom


####################################################################################################################################
#####       Optimized (single fraction) Humlicek region I rational approximation for n=16 and delta=1.35                       #####
####################################################################################################################################

AA = np.array([  +46236.3358828121,
                 -147726.58393079657j,
                 -206562.80451354137,
                 281369.1590631087j,
                 183092.74968253175,
                 -184787.96830696272j,
                 -66155.39578477248,
                 57778.05827983565j,
                 11682.770904216826,
                 -9442.402767960672j,
                 -1052.8438624933142,
                 814.0996198624186j,
                 45.94499030751872,
                 -34.59751573708725j,
                 -0.7616559377907136,
                 0.5641895835476449j])  # identical to 1/sqrt(pi) for 12 digits, except for last four digits

# NOTE:  BB=bb --- the denominator coefficients are identical to those for delta=1.31183, i.e. independent of delta!!!

numPoly16a = np.poly1d(np.flipud(AA))

def zpf16P (z):
	""" Humlicek (1979) complex probability function for n=16 and delta=1.35.
	    Optimized rational approximation using numpy.poly1d  (applicable for all z).

	    Maximum relative error to wofz in 0<x<25 and 1e-8<y<1e2:  1.992e-4 and 2.60e-6 for K and L
	    """
	Z = z + 1.35j
	return numPoly16a(Z)/denPoly8e(Z*Z)


####################################################################################################################################
#####  Combinations of optimized Humlicek (1979) rational approximation (n=16, delta=1.35) with his 1982 asymptotic approx.    #####
####################################################################################################################################

def hum1zpf16w (z, s=15.0):
	""" Complex error function w(z)=w(x+iy) combining Humlicek's rational approximations:

	    |x|+y>s:   Humlicek (JQSRT, 1982) rational approximation for region I;
	    else:      Humlicek (JQSRT, 1979) rational approximation with n=16 and delta=y0=1.35

	    Version using np.where, useful for contour plots.

	    Maximum relative error to wofz in 0<x<25 and 1e-8<y<1e2:  1.462e-4 and 7.73e-5 for K and L
	                                                (1e-6<y<1e2:  1.393e-5 for K with s=25)
	"""

	return np.where(abs(z.real)+z.imag>s, z * iRecSqrtPi / (z*z-0.5), zpf16P(z))


def hum2zpf16w (z, s=6.0):
	""" Complex error function w(z)=w(x+iy) combining Humlicek's rational approximations:

	    |x|+y>s:   Humlicek (JQSRT, 1982) rational approximation for region II;
	    else:      Humlicek (JQSRT, 1979) rational approximation with n=16 and delta=y0=1.35

	    Version using np.where, useful for contour plots.

	    Maximum relative error to wofz in 0<x<25 and 1e-8<y<1e2:  7.310e-5 and 1.74e-5 for K and L
	    """

	zz = z*z
	return np.where(abs(z.real)+z.imag>s,
	                1j* (z * (zz*recSqrtPi-1.410474))/ (0.75 + zz*(zz-3.0)),
	                zpf16P(z))


####################################################################################################################################
#####  Humlicek (1979, 1982) combinations as above, but with real arguments x,y and Horner scheme for zpf16                    #####
#####  !! these versions are more useful for line-by-line 'number crunching' (each line has its own xGrid, but a common y) !!  #####
####################################################################################################################################

def hum1zpf16m (x, y, s=15.0):
	""" Complex error function w(z)=w(x+iy) combining Humlicek's rational approximations:

	    |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
	    else:      Humlicek (JQSRT, 1979) rational approximation with n=16 and delta=y0=1.35

	    Version using a mask and np.place;  two real arguments x,y.  """

	# For safety only. Probably slows down everything. Comment it if you always have arrays (or use assertion?).
	# if isinstance(x,(int,float)):  x = np.array([x])

	t = y-1j*x
	w = t * recSqrtPi / (0.5 + t*t)  # Humlicek (1982) approx 1 for s>15

	if y<s:
		mask  = abs(x)+y<s                      # returns true for interior points
		Z     = x[np.where(mask)]+ 1j*(y+1.35)  # returns small complex array covering only the interior region
		ZZ    = Z*Z
		numer = (((((((((((((((AA[15]*Z+AA[14])*Z+AA[13])*Z+AA[12])*Z+AA[11])*Z+AA[10])*Z+AA[9])*Z+AA[8])*Z+AA[7])*Z+AA[6])*Z+AA[5])*Z+AA[4])*Z+AA[3])*Z+AA[2])*Z+AA[1])*Z+AA[0])
		denom = (((((((ZZ+bb[14])*ZZ+bb[12])*ZZ+bb[10])*ZZ+bb[8])*ZZ+bb[6])*ZZ+bb[4])*ZZ+bb[2])*ZZ+bb[0]
		np.place(w, mask, numer/denom)
	return w


def hum2zpf16m (x, y, s=10.0):
	""" Complex error function w(z)=w(x+iy) combining Humlicek's rational approximations:

	    |x|+y>10:  Humlicek (JQSRT, 1982) rational approximation for region II;
	    else:      Humlicek (JQSRT, 1979) rational approximation with n=16 and delta=y0=1.35

	    Version using a mask and np.place;  two real arguments x,y. """

	# For safety only. Probably slows down everything. Comment it if you always have arrays (or use assertion?).
	# if isinstance(x,(int,float)):  x = np.array([x])

	z  = x+1j*y
	zz = z*z
	w  = 1j* (z * (zz*recSqrtPi-1.410474))/ (0.75 + zz*(zz-3.0))

	if y<s:
		mask  = abs(x)+y<s                      # returns true for interior points
		Z     = x[np.where(mask)]+ 1j*(y+1.35)  # returns small complex array covering only the interior region
		ZZ    = Z*Z
		numer = (((((((((((((((AA[15]*Z+AA[14])*Z+AA[13])*Z+AA[12])*Z+AA[11])*Z+AA[10])*Z+AA[9])*Z+AA[8])*Z+AA[7])*Z+AA[6])*Z+AA[5])*Z+AA[4])*Z+AA[3])*Z+AA[2])*Z+AA[1])*Z+AA[0])
		denom = (((((((ZZ+bb[14])*ZZ+bb[12])*ZZ+bb[10])*ZZ+bb[8])*ZZ+bb[6])*ZZ+bb[4])*ZZ+bb[2])*ZZ+bb[0]
		np.place(w, mask, numer/denom)
	return w


####################################################################################################################################
#####         Time benchmarks       (November 2022,  Intel(R) Core(TM) i5-9600 CPU 3.10GHz,  Suse 15.3)                        #####
#####                                                                                                                          #####
#####         (see also Tab. 2 and Fig. 5 in the MNRAS 2018 paper)                                                             #####
####################################################################################################################################

# In [1]: from cpfX import *

# In [2]: x=np.linspace(0.,25.,251)
#    ...: y=np.logspace(-8,2,1001)
#    ...: # kind of 'matrices' required for the contour plots (note Python is case sensitive!)
#    ...: X, Y = np.meshgrid(x, y)
#    ...: Z = X +1j*Y

# In [3]: %timeit zpf16l(Z)
# 15.2 ms plusminus 0.0184 ms per loop (mean plusminus std. dev. of 7 runs, 100 loops each)

# In [4]: %timeit zpf16h(Z)
# 11.3 ms plusminus 0.0232 ms per loop (mean plusminus std. dev. of 7 runs, 100 loops each)

# In [5]: %timeit zpf16p(Z)
# 16.4 ms plusminus 0.029 ms per loop (mean plusminus std. dev. of 7 runs, 100 loops each)

# In [6]: x=np.linspace(0.,25.,2501)
#    ...: y=np.logspace(-8,2,1001)
#    ...: # kind of 'matrices' required for the contour plots (note Python is case sensitive!)
#    ...: X, Y = np.meshgrid(x, y)
#    ...: Z = X +1j*Y

# In [7]: %timeit zpf16l(Z)
# 231 ms plusminus 1.48 ms per loop (mean plusminus std. dev. of 7 runs, 1 loop each)

# In [8]: %timeit zpf16h(Z)
# 172 ms plusminus 0.347 ms per loop (mean plusminus std. dev. of 7 runs, 10 loops each)

# In [9]: %timeit zpf16p(Z)
# 244 ms plusminus 1.24 ms per loop (mean plusminus std. dev. of 7 runs, 1 loop each)
