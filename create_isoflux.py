
import numpy
import matplotlib.pyplot
def create_isoflux(npts, r0, z0, a, kappa, delta,  zeta = 0, upnull = False, lonull = False, kappa_L = None, delta_L = None, zeta_IU = None, zeta_IL = None, zeta_OL = None):

    if zeta_IU is None: zeta_IU = zeta
    if zeta_IL is None: zeta_IL = zeta
    if zeta_OL is None: zeta_OL = zeta
    if kappa_L is None: kappa_L = kappa
    if delta_L is None: delta_L = delta

    ang1 = numpy.linspace(0,numpy.pi/2, npts, endpoint = False)
    ang2 = numpy.linspace(numpy.pi/2,2*numpy.pi/2, npts, endpoint = False)
    ang3 = numpy.linspace(2*numpy.pi/2,3*numpy.pi/2, npts, endpoint = False)
    ang4 = numpy.linspace(3*numpy.pi/2,4*numpy.pi/2, npts, endpoint = False)

    zetas = numpy.array([zeta, zeta_IU, zeta_IL, zeta_OL])
    zetas[zetas < -1.0 * 1.0 / numpy.sqrt(2.0)] = -1.0 * 1.0 / numpy.sqrt(2.0)

    n1 = -numpy.log(2.0) / numpy.log(zetas[0] * (1.0 - 1.0 / numpy.sqrt(2.0)) + 1.0 / numpy.sqrt(2.0))
    r1 = a * (1.0 / (a/r0) - delta) + a * (1.0 + delta) * numpy.cos(ang1) ** (2.0 / n1)
    z1 = z0 + a * kappa * numpy.sin(ang1) ** (2.0 / n1)
    n2 = -numpy.log(2.0) / numpy.log(zetas[1] * (1.0 - 1.0 / numpy.sqrt(2.0)) + 1.0 / numpy.sqrt(2.0))
    r2 = a * (1.0 / (a/r0) - delta) - a * (1.0 - delta) * abs(numpy.cos(ang2)) ** (2.0 / n2)
    z2 = z0 + a * kappa * numpy.sin(ang2) ** (2.0 / n2)
    n3 = -numpy.log(2.0) / numpy.log(zetas[2] * (1.0 - 1.0 / numpy.sqrt(2.0)) + 1.0 / numpy.sqrt(2.0))
    r3 = a * (1.0 / (a/r0) - delta_L) - a * (1.0 - delta_L) * abs(numpy.cos(ang3)) ** (2.0 / n3)
    z3 = z0 - a * kappa_L * abs(numpy.sin(ang3)) ** (2.0 / n3)
    n4 = -numpy.log(2.0) / numpy.log(zetas[3] * (1.0 - 1.0 / numpy.sqrt(2.0)) + 1.0 / numpy.sqrt(2.0))
    r4 = a * (1.0 / (a/r0) - delta_L) + a * (1.0 + delta_L) * abs(numpy.cos(ang4)) ** (2.0 / n4)
    z4 = z0 - a * kappa_L * abs(numpy.sin(ang4)) ** (2.0 / n4)

    a_s = a*numpy.array([1+delta, 1-delta, 1-delta_L, 1+delta_L])
    b_s = a*numpy.array([kappa, kappa, kappa_L, kappa_L])
    g_s = numpy.array([delta - 1 if delta >= 0 else -1.0 / (1.0 + delta), -1.0 / (1.0 - delta) if delta >= 0.0 else -1.0 * (1.0 + delta), -1.0 / (1.0 - delta_L) if delta_L >= 0.0 else -1.0 * (1.0 + delta_L), delta_L - 1 if delta_L >= 0 else -1.0 / (1.0 + delta_L)])
    h_s = 1 - (1-zetas)*(1-1/numpy.sqrt(2))

    f_s = numpy.zeros(4)
    n_s = numpy.ones(4)

    for i in range(4):
        if ((i==0 or i==1) and upnull) or ((i==2 or i==3) and lonull):
            f = numpy.linspace(0.01, 0.99, 100)
            n = numpy.linspace(1.04, 5, 100)
            y1 = numpy.zeros((100, 100))
            y2 = numpy.zeros((100, 100))
            for j in range(y2.shape[0]):
                for k in range(y1.shape[1]):
                    y1[k, j] = (f[j] + h_s[i] * (1.0 - f[j])) ** n[k] + (1.0 - f[j] ** n[k]) * h_s[i] ** n[k] - 1.0
                    y2[k, j] = f[j] ** (n[k] - 1.0) * (f[j] * (g_s[i] + b_s[i] / a_s[i]) - b_s[i]/a_s[i]) - g_s[i]
            fig, ax = matplotlib.pyplot.subplots()
            xy1 = ax.contour(f, n, y1, levels = [0])
            xy1 = numpy.array(xy1.get_paths()[0].vertices)[numpy.argsort(numpy.array(xy1.get_paths()[0].vertices[:, 0]))]
            xy2 = ax.contour(f, n, y2, levels = [0])
            xy2 = numpy.array(xy2.get_paths()[0].vertices)[numpy.argsort(numpy.array(xy2.get_paths()[0].vertices[:, 0]))]
            matplotlib.pyplot.close()
            y1sol = numpy.interp(f, xy1[:, 0], xy1[:, 1])
            y2sol = numpy.interp(f, xy2[:, 0], xy2[:, 1])
            maxdiff = max(y1sol - y2sol)
            mindiff = min(y1sol - y2sol)
            if maxdiff / mindiff < 0.0:
                y12sol = min(abs(y1sol - y2sol))
                imin = numpy.argmin(abs(y1sol - y2sol))
                f_s[i] = f[imin]
                n_s[i] = y1sol[imin]
            else:
                if maxdiff > 0:
                    y1new = (f[94] + h_s[i] * (1.0 - f[94])) ** y2sol[94] + (1.0 - f[94] ** y2sol[94]) * h_s[i] ** y2sol[94] - 1.0
                    y2new = f[94] ** (y2sol[94] - 1.0) * (f[94] * (g_s[i] + b_s[i] / a_s[i]) - b_s[i] / a_s[i]) - g_s[i]
                    while y1new > y2new:
                        h_s[i] = h_s[i] - 0.01
                        y1new = (f[94] + h_s[i] * (1.0 - f[94])) ** y2sol[94] + (1.0 - f[94] ** y2sol[94]) * h_s[i] ** y2sol[94] - 1.0
                    f_s[i] = f[94]
                    n_s[i] = y2sol[94]
                else:
                    y1new = (f[4] + h_s[i] * (1.0 - f[4])) ** y2sol[4] + (1.0 - f[4] ** y2sol[4]) * h_s[i] ** y2sol[4] - 1.0
                    y2new = f[4] ** (y2sol[4] - 1.0) * (f[4] * (g_s[i] + b_s[i] / a_s[i]) - b_s[i] / a_s[i]) - g_s[i]
                    while y1new < y2new:
                        h_s[i] = h_s[i] + 0.01
                        y1new = (f[4] + h_s[i] * (1.0 - f[4])) ** y2sol[4] + (1.0 - f[4] ** y2sol[4]) * h_s[i] ** y2sol[4] - 1.0
                    f_s[i] = f[4]
                    n_s[i] = y2sol[4]
    zetas_new = 1-(1-h_s)/(1-1/numpy.sqrt(2))
    alphas = a_s/(1-f_s)
    betas = b_s/((1-f_s**n_s)**(1/n_s))

    if upnull:
        z1 = z0 + betas[0] * (1.0 - ((r1 - a * (1.0 / (a/r0) + 1.0)) / alphas[0] + 1.0) ** n_s[0]) ** (1.0 / n_s[0])
        z2 = z0 + betas[1] * (1.0 - (1.0 + (a * (1.0 / (a/r0) - 1.0) - r2) / alphas[1]) ** n_s[1]) ** (1.0 / n_s[1])
    if lonull:
        z3 = z0-betas[2] * (1.0 - (1.0 + (a * (1.0 / (a/r0) - 1.0) - r3) / alphas[2]) ** n_s[2]) ** (1.0 / n_s[2])
        z4 = z0 -betas[3] * (1.0 - ((r4 - a * (1.0 / (a/r0) + 1.0)) / alphas[3] + 1.0) ** n_s[3]) ** (1.0 / n_s[3])

    r = numpy.hstack((r1, r2, r3, r4))
    z = numpy.hstack((z1, z2, z3, z4))
    return numpy.column_stack([r, z]), zetas_new

