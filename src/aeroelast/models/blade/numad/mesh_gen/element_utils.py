# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 12:44:35 2023

@author: evans
"""

import numpy as np


def cross_prod(v1, v2):
    cp = np.zeros(3, dtype=float)
    cp[0] = v1[1] * v2[2] - v1[2] * v2[1]
    cp[1] = v1[2] * v2[0] - v1[0] * v2[2]
    cp[2] = v1[0] * v2[1] - v1[1] * v2[0]
    return cp


def rotate_orient(oriMat, rotVec):
    mag = np.linalg.norm(rotVec)
    if mag > 0.0001:
        a1 = np.zeros((3, 3), dtype=float)
        sth = np.sin(mag)
        cth = np.cos(mag)

        unitRot = (1.0 / mag) * rotVec
        a1[0] = unitRot

        i1 = 0
        if abs(unitRot[1]) < abs(unitRot[0]):
            i1 = 1
        if abs(unitRot[2]) < abs(unitRot[i1]):
            i1 = 2

        a1[1, i1] = np.sqrt(1.0 - unitRot[i1] * unitRot[i1])

        for i2 in range(0, 3):
            if i2 != i1:
                a1[1, i2] = -a1[0, i1] * a1[0, i2] / a1[1, i1]

        a1[2] = cross_prod(a1[0], a1[1])

        a2 = np.array([[1.0, 0.0, 0.0], [0.0, cth, -sth], [0.0, sth, cth]])

        a3 = np.matmul(a2, a1)
        a3T = np.transpose(a3)
        a4 = np.matmul(a3T, a1)

        return np.matmul(oriMat, a4)
    else:
        return oriMat


def correct_orient(aprxOri, elCrd, elType):
    if elType == "shell4":
        v0 = elCrd[:, 2] - elCrd[:, 0]
        v1 = elCrd[:, 3] - elCrd[:, 1]
    elif elType == "shell3":
        v0 = elCrd[:, 1] - elCrd[:, 0]
        v1 = elCrd[:, 2] - elCrd[:, 0]
    else:
        pstr = "Warning: element type " + elType + " not supported for correctOrient"
        print(pstr)
        return aprxOri
    v2 = cross_prod(v0, v1)
    mag = np.linalg.norm(v2)
    v2 = (1.0 / mag) * v2
    a2 = aprxOri[2]
    dp = np.dot(v2, a2)
    if dp < 0.0:
        v2 = -1.0 * v2
    cp = cross_prod(a2, v2)
    mag = np.linalg.norm(cp)
    theta = np.arcsin(mag)
    rotVec = (theta / mag) * cp
    return rotate_orient(aprxOri, rotVec)
