"""V2-only numerical tools; parsed float inputs keep their binary values.

Exact rational chart predicates are separate from extended Cartesian arithmetic.
The latter requires NumPy longdouble to have more precision than binary64 and
fails closed when final vertices/Jacobians cannot represent a finite band.
No epsilon identifies source events or boundaries.
"""
from fractions import Fraction

import numpy as np


DTYPE = np.longdouble


def extended(value):
    """Convert a rational without the implicit Fraction -> binary64 detour."""
    if isinstance(value, Fraction):
        return DTYPE(str(value.numerator))/DTYPE(str(value.denominator))
    return DTYPE(value)


def extended_array(values):
    array = np.asarray(values, dtype=object)
    return np.array([extended(v) for v in array.flat], dtype=DTYPE).reshape(array.shape)


def rational(value):
    """Retain the exact binary value of floats, including longdouble scalars."""
    if isinstance(value, Fraction):
        return value
    return Fraction(*value.as_integer_ratio()) if hasattr(value, 'as_integer_ratio') else Fraction(value)


def linspace(a, b, count):
    """Exact chart subdivision, including both endpoints."""
    a, b = rational(a), rational(b)
    return [a+(b-a)*i/(count-1) for i in range(count)]


class Pchip:
    """Fritsch-Butland PCHIP in longdouble, axis zero, value and derivative.

    Same shape-preserving formula as SciPy PchipInterpolator; SciPy coerces its
    coefficients/queries to binary64, so cannot evaluate the extended v2 path.
    This is not a new interpolation policy or a higher-order loft.
    """
    def __init__(self, x, y, axis=0):
        if axis != 0:
            raise ValueError('Extended windIO PCHIP supports axis=0 only')
        if np.finfo(DTYPE).nmant <= np.finfo(float).nmant:
            raise ValueError('windIO extended loft requires longdouble wider than float64')
        self.x = extended_array(x)
        y = extended_array(y)
        h = np.diff(self.x)
        if len(h) < 1 or np.any(h <= 0) or y.shape[0] != len(self.x):
            raise ValueError('Extended windIO PCHIP requires increasing knots')
        h = h.reshape((-1,)+(1,)*(y.ndim-1))
        slopes = np.diff(y, axis=0)/h
        d = np.zeros_like(y)
        if len(y) == 2:
            d[:] = slopes[0]
        else:
            before, after = slopes[:-1], slopes[1:]
            same = (before != 0) & (after != 0) & (np.sign(before) == np.sign(after))
            w1, w2 = 2*h[1:]+h[:-1], h[1:]+2*h[:-1]
            safe_before, safe_after = np.where(same, before, 1), np.where(same, after, 1)
            d[1:-1] = np.where(same, (w1+w2)/(w1/safe_before+w2/safe_after), 0)

            def endpoint(h0, h1, m0, m1):
                result = ((2*h0+h1)*m0-h0*m1)/(h0+h1)
                result = np.where(np.sign(result) != np.sign(m0), 0, result)
                return np.where((np.sign(m0) != np.sign(m1)) & (abs(result) > 3*abs(m0)),
                                3*m0, result)
            d[0] = endpoint(h[0], h[1], slopes[0], slopes[1])
            d[-1] = endpoint(h[-1], h[-2], slopes[-1], slopes[-2])
        self.y, self.d = y, d
        self.c2 = (3*slopes-2*d[:-1]-d[1:])/h
        self.c3 = (d[:-1]+d[1:]-2*slopes)/(h*h)

    def __call__(self, query, derivative=0):
        q = extended_array(query)
        indices = np.clip(np.searchsorted(self.x, q, side='right')-1, 0, len(self.x)-2)
        t = (q-self.x[indices]).reshape(q.shape+(1,)*(self.y.ndim-1))
        endpoint = (q == self.x[indices+1]).reshape(q.shape+(1,)*(self.y.ndim-1))
        if derivative == 1:
            result = (3*self.c3[indices]*t+2*self.c2[indices])*t+self.d[indices]
            return np.where(endpoint, self.d[indices+1], result)
        result = ((self.c3[indices]*t+self.c2[indices])*t+self.d[indices])*t+self.y[indices]
        return np.where(endpoint, self.y[indices+1], result)

    def derivative(self):
        return lambda query: self(query, derivative=1)

    def directional(self, query, dy):
        """Derivative with respect to ordinates moving along ``dy``.

        Used for the loft's analytic s derivative: interpolation in rthick
        depends nonlinearly on the master ordinates at s. Differentiating
        those slope formulas avoids finite differences and noisy minimizers.
        """
        h = np.diff(self.x).reshape((-1,)+(1,)*(self.y.ndim-1))
        m, dm = np.diff(self.y, axis=0)/h, np.diff(extended_array(dy), axis=0)/h
        dd = np.zeros_like(self.d)
        if len(self.y) == 2:
            dd[:] = dm[0]
        else:
            same = (m[:-1] != 0) & (m[1:] != 0) & (np.sign(m[:-1]) == np.sign(m[1:]))
            before, after = np.where(same, m[:-1], 1), np.where(same, m[1:], 1)
            w1, w2 = 2*h[1:]+h[:-1], h[1:]+2*h[:-1]
            dd[1:-1] = np.where(same, self.d[1:-1]**2/(w1+w2)*(
                w1*dm[:-1]/before**2+w2*dm[1:]/after**2), 0)

            def endpoint(h0, h1, m0, m1, dm0, dm1):
                raw = ((2*h0+h1)*m0-h0*m1)/(h0+h1)
                derivative = ((2*h0+h1)*dm0-h0*dm1)/(h0+h1)
                clipped = (np.sign(m0) != np.sign(m1)) & (abs(raw) > 3*abs(m0))
                return np.where(np.sign(raw) != np.sign(m0), 0,
                                np.where(clipped, 3*dm0, derivative))
            dd[0] = endpoint(h[0], h[1], m[0], m[1], dm[0], dm[1])
            dd[-1] = endpoint(h[-1], h[-2], m[-1], m[-2], dm[-1], dm[-2])
        tangent = object.__new__(Pchip)
        tangent.x, tangent.y, tangent.d = self.x, extended_array(dy), dd
        tangent.c2 = (3*dm-2*dd[:-1]-dd[1:])/h
        tangent.c3 = (dd[:-1]+dd[1:]-2*dm)/(h*h)
        return tangent(query)


def minimum_location(derivative):
    """Locate the single-LE loft's stationary point in extended arithmetic.

    Bisection of the analytic derivative converges in position, unlike an
    objective-value minimizer whose flat minimum loses about half the bits.
    This selects samples only; it neither moves nor redefines the source loft.
    """
    # TE extrema need not bracket the interior LE minimum. Isolate the
    # negative-to-positive transition, excluding positive-to-negative maxima.
    # This is a bracket search only: no source coordinate is moved or merged.
    samples = np.linspace(DTYPE(0), DTYPE(1), 65)
    slopes = np.asarray(derivative(samples), dtype=DTYPE)
    if slopes.shape != samples.shape or not np.isfinite(slopes).all():
        raise ValueError('windIO surface leading-edge landmark is not bracketed')
    nonzero = np.flatnonzero(slopes != 0)
    brackets = [(samples[i], samples[j])
                for i, j in zip(nonzero[:-1], nonzero[1:], strict=True)
                if slopes[i] < 0 < slopes[j]]
    if len(brackets) != 1:
        raise ValueError('windIO surface leading-edge landmark is not bracketed')
    a, b = brackets[0]
    for _ in range(np.finfo(DTYPE).nmant+4):
        middle = a+(b-a)/2
        if not a < middle < b:
            break
        slope = derivative(middle)
        if slope == 0:
            return middle
        if slope < 0:
            a = middle
        else:
            b = middle
    return a+(b-a)/2
