# -*- coding: utf-8 -*-
"""
Created on Sun Nov 08 01:03:05 2015

@author: Kevin Multani (kvmu)

Description: TO DO

Original code by: Yuya Horita
"""
import numpy as np

def coordTransform(data, coordinate_system):
    if coordinate_system == 'xyz':
        return data

    elif coordinate_system == 'xzy':
        data[:, 1], data[:, 2] = -data[:, 2], data[:, 1]
        data[:, 4], data[:, 5] = -data[:, 5], data[:, 4]
        return data

    elif coordinate_system == 'yxz':
        data[:, 0], data[:, 1] = data[:, 1], -data[:, 0]
        data[:, 4], data[:, 5] = data[:, 5], -data[:, 4]
        return data

    else:
        print 'invalid coordinate system'
        raise ValueError

def get_x_y_z_data(data):
    return (data[:, 0], data[:, 1], data[:, 2],
            data[:, 3], data[:, 4], data[:, 5])


def to_r_theta_z_position(x, y, z):
    r = np.sqrt(np.square(x) + np.square(y))
    theta = np.rad2deg(np.arctan(np.divide(y, x)))
    z = z
    return r, theta, z


def to_r_theta_z_field(theta, Bx, By, Bz):
    rad = np.deg2rad(theta)
    Br = Bx * np.cos(rad) + By * np.sin(rad)
    Btheta = -Bx * np.sin(rad) + By * np.cos(rad)
    Bz = Bz
    return Br, Btheta, Bz


def getData(filename, coordinate_system, symmetry_number):
    data_matrix = np.genfromtxt(filename, skip_header=8)
    data_matrix = coordTransform(data_matrix, coordinate_system)
    x, y, z, Bx, By, Bz = get_x_y_z_data(data_matrix)
    r, theta, z = to_r_theta_z_position(x, y, z)
    Br, Btheta, Bz = to_r_theta_z_field(theta, Bx, By, Bz)

    # Assert mid-plane symettery
    Br[np.abs(z) < 1e-9] = 0
    Btheta[np.abs(z) < 1e-9] = 0
    return (np.round(r, decimals=5),
            np.round(theta, decimals=5),
            np.round(z, decimals=5),
            Br, Btheta, Bz)

def getSize(r, theta, z):
    return (np.where(r != r[0])[0][0],
            np.where(theta != theta[0])[0][0],
            np.where(z != z[0])[0][0])

def getAveFields(r, theta, z, Br, Btheta, Bz,
                 r_size, theta_size, z_size,
                 symmetry_number):
    return








