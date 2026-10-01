from scipy.interpolate import PchipInterpolator, interp1d


def interpolator_wrap(x, v, xq, method="linear", axis=0):
    """This function is designed to emulate the arg structure and output
    of matlabs interp1d function.

    Parameters
    ----------
        x : array
        v : array
        xq : array
        method : str
            Defaults to 'linear'.
        axis : int
            Defaults to 0.
        extrapolation :bool
            Defaults to None.

    Returns:
        array :
    """
    if method == "linear":
        interpolator = interp1d(x, v, "linear", axis, bounds_error=False, fill_value="extrapolate")
        vq = interpolator(xq)
    elif method == "pchip":
        interpolator = PchipInterpolator(x, v, axis, extrapolate=True)
        vq = interpolator(xq)
    elif method == "spline":
        interpolator = interp1d(x, v, "cubic", axis, bounds_error=False, fill_value="extrapolate")
        vq = interpolator(xq)
    return vq
    # if method == 'pp':
    #     pass
    # if method == 'v5cubic':
    #     raise Exception("Method error for interpolator_wrap. 'v5cubic' not implemented")
    # if method == 'makima':
    #     raise Exception("Method error for interpolator_wrap. 'makima' not implemented")
    # if method == 'nearest':
    #     raise Exception("Method error for interpolator_wrap. 'nearest' not implemented")
    # if method == 'next':
    #     raise Exception("Method error for interpolator_wrap. 'next' not implemented")
    # if method == 'previous':
    #     raise Exception("Method error for interpolator_wrap. 'previous' not implemented")
