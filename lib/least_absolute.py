import numpy as np


def least_absolute(x, y):  # Least Absolute
    """
    @param x: a value list
    @param y: a value list
    @return: final = [a, b] // y = a*x + b
    """
    xo, yo = np.median(x), np.median(y)
    a, b = [], []

    nx, ny = len(x), len(y)
    n = nx
    if nx != ny:
        print'\nSorry! The data is not twinning, please check your data ^o^ !\n'

    for i in range(n):
        dealtx, dealty = x[i] - xo, y[i] - yo
        if dealtx != 0:
            k = dealty / dealtx
            a.append(k)
        else:
            a.append(10000)

    for i in range(n):
        bb = y[i] - a[i] * x[i]
        b.append(bb)

    error = []
    err = []
    for i in range(n):
        for k in range(n):
            ee = abs((a[i] * x[k] + b[i]) - y[k])
            error.append(ee)
        temp = sum(error)

        err.append(temp)
        error = []
    minimum = min(err)

    place = []
    for i, val in enumerate(err):
        if val == minimum:
            place.append(i)

    final = [a[place[0]], b[place[0]]]

    kspace, bspace = [], []
    if n > 100:
        nk = 100
        nstep = 0.00001
        nprint = 5000
    else:
        nk = 1000
        nstep = 0.001
        nprint = 100000
    nb = nk + 1

    for i in range(nk):
        j = nk - i
        tk = final[0] - nstep * j
        kspace.append(tk)

    for i in range(nb):
        ttk = final[0] + nstep * i
        kspace.append(ttk)

    for i in range(nk):
        j = nk - i
        tb = final[1] - nstep * j
        bspace.append(tb)

    for i in range(nb):
        ttb = final[1] + nstep * i
        bspace.append(ttb)

    lk, lb = len(kspace), len(bspace)
    ff, ferr = [], []
    count = 0
    for i in range(lk):
        for j in range(lb):
            for k in range(n):
                temp6 = kspace[i] * x[k] + bspace[j]
                temp7 = abs(temp6 - y[k])
                ff.append(temp7)
            ferr.append(sum(ff))
            ff = []
            count += 1

    minimum = min(ferr)

    pla = []
    for i, val in enumerate(ferr):
        if val == minimum:
            pla.append(i)

    site = divmod(pla[0], nk * 2 + 1)
    refinal = [kspace[site[0]], bspace[site[1] - 1]]

    return refinal
