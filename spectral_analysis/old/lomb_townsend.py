import numpy as np
 
def lomb(x,y,freqs):
	# Check input sizes
    if x.shape[0] != y.shape[0]:
	    raise ValueError("Input arrays do not have the same size.")

    # Create empty array for output periodogram
    pgram = np.empty(freqs.shape[0], dtype=np.float64)
    ry = np.empty(freqs.shape[0], dtype=np.float64)
    iy = np.empty(freqs.shape[0], dtype=np.float64)

    vari = np.var(y)

    # Local variables
	#cdef Py_ssize_t i, j
	#cdef double c, s, xc, xs, cc, ss, cs
	#cdef double tau, c_tau, s_tau, c_tau2, s_tau2, cs_tau

    for i in range(freqs.shape[0]):

        xc = 0.
        xs = 0.
        cc = 0.
        ss = 0.
        cs = 0.

        for j in range(x.shape[0]):

            c = np.cos(freqs[i] * x[j])
            s = np.sin(freqs[i] * x[j])
            
            xc += y[j] * c
            xs += y[j] * s
            cc += c * c
            ss += s * s
            cs += c * s

        tau = np.arctan(2 * cs / (cc - ss)) / (2 * freqs[i])
        c_tau = np.cos(freqs[i] * tau)
        s_tau = np.sin(freqs[i] * tau)
        c_tau2 = c_tau * c_tau
        s_tau2 = s_tau * s_tau
        cs_tau = 2 * c_tau * s_tau

        ry[i] = ((c_tau * xc + s_tau * xs)**2 / (c_tau2 * cc + cs_tau * cs + s_tau2 * ss))
        iy[i] = ((c_tau * xs - s_tau * xc)**2 / (c_tau2 * ss - cs_tau * cs + s_tau2 * cc))

        PHLS=np.arctan2(iy,ry)

        pgram[i] = 0.5 * (ry[i] + iy[i])/vari

    dimi=2.*(len(x)+1.)
    fac=np.sqrt(vari*dimi/2.)
    amp=fac*np.sqrt(pgram)    

    corr_amp = amp/len(x)
    corr_amp = corr_amp*2.

    return pgram,corr_amp
