from math import log, exp
import numpy as np
from numpy.ma import log2
from scipy.misc import comb


def cal_mavalue(pk, mk1_name, mk2_name):
    """
    calculate M&A value of known pk with 2 reads data
    """
    if len(pk.rds_density) >= 2 and mk1_name in pk.rds_density and mk2_name in pk.rds_density:
        density1 = pk.rds_density[mk1_name]
        density2 = pk.rds_density[mk2_name]
        mvalue = log2(density1) - log2(density2)
        avalue = (log2(density1) + log2(density2)) / 2
        pk.another_info.update({'mvalue': mvalue})
        pk.another_info.update({'avalue': avalue})
        return mvalue, avalue


def cal_mapvalue_rescaled(pk, mk1_name, mk2_name, ma_fit):
    """
    calculate M&A&P value of known pk with 2 reads data and fit parameters
    ma_fit: R2 = ma_fit[0] * R1 + ma_fit[1]
    """
    density1 = pk.rds_density[mk1_name]
    density2 = pk.rds_density[mk2_name]
    log2_density1_re = (2 - ma_fit[1]) * log2(density1) / (2 + ma_fit[1]) - 2 * ma_fit[0] / (2 + ma_fit[1])
    mvalue_re = log2_density1_re - log2(density2)
    avalue_re = (log2_density1_re + log2(density2)) / 2

    density1_norm = 2 ** log2_density1_re
    density2_norm = 2 ** log2(density2)
    pvalue = np.ones(pk.pk_num)
    for i in xrange(pk.pk_num):
        pvalue[i] = __digit_exprs_p_norm(density1_norm[i], density2_norm[i])
    pk.another_info.update({'MAnorm_mvalue': mvalue_re})
    pk.another_info.update({'MAnorm_avalue': avalue_re})
    pk.another_info.update({'MAnorm_pvalue': pvalue})
    return mvalue_re, avalue_re, pvalue


def __digit_exprs_p_norm(x, y):
    xx = round(x)
    if xx == 0:
        xx = 1
    yy = int(round(y))
    if xx + yy < 20.0:  # if x + y small
        p1 = round(comb(xx + yy, xx)) * 2 ** - (xx + yy + 1.0)
        p2 = round(comb(xx + yy, yy)) * 2 ** - (xx + yy + 1.0)
        return max(p1, p2)
    else:  # if x + y large, use the approximate equations
        log_p = (xx + yy) * log(xx + yy) - xx * log(xx) - yy * log(yy) - (xx + yy + 1.0) * log(2.0)
        if log_p < -500:
            log_p = -500
        p = exp(log_p)
        return p


if __name__ == '__main__':
    v = __digit_exprs_p_norm(64, 1492)
    pass