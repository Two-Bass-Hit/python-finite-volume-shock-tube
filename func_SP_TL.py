# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:48:12 2019

@author: j_karlovsky
"""
import numpy as np


def SP_TL(area_1, area_2, mag_spec_1, mag_spec_2):

    numerator = area_1 * pow(mag_spec_1, 2)  #raise mag_spec_1 vector to the power 2 (square it)
    denom = area_2 * pow(mag_spec_2, 2)
    TL = 10 * np.log10(numerator/denom)
    return TL