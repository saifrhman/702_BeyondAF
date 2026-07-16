# -*- coding = utf-8 -*-
# @File: math_base.PY
from typing import Tuple, Union

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.floating]


def get_distance(data_list):
    """

    :param data_list: A vector
    :return: Value indicates the length of vector
    """
    data_list = [i**2 for i in data_list]
    return sum(data_list) ** 0.5


def get_strength(data_list):
    """Compute the strength of the triangle sigma=(p-a)(p-b)(p-c)/p^2, where p=(a+b+c)/2.

    :param data_list: List containing length of sides of a polygon
    :return: Strength calculated with the formula above
    """
    p = sum(data_list) * 0.5
    return np.prod([p - i for i in data_list]) / (p**2)


def get_norm_strength(data_list) -> float:
    """Compute the strength of the triangle sigma=(p-a)(p-b)(p-c)/p^2, where p=(a+b+c)/2.

    :param data_list: List containing length of sides of a polygon
    :return: Strength calculated with the formula above
    """
    p = sum(data_list) * 0.5
    return 27 * np.prod([p - i for i in data_list]) / (p**3)


def dot_product(v1: FloatArray, v2: FloatArray) -> float:
    """Compute dot product of 2 vectors

    :param v1: 1-D array vector
    :param v2: 1-D array vector
    :return: Dot product of input vectors
    """
    return np.dot(v1, v2)


def cross_product(v1: FloatArray, v2: FloatArray) -> float:
    """Compute cross product of 2 vectors

    :param v1: 1-D array vector
    :param v2: 1-D array vector
    :return: Cross product of input vectors
    """
    return np.cross(v1, v2)


def vector_norm(vector: FloatArray) -> float:
    """Normalise the input vector

    :param vector: 1-D array vector
    :return: Length of the input vector
    """
    return np.linalg.norm(vector)


def vector_round(n: FloatArray, decimal: int = 3) -> FloatArray:
    """Round the input vetor to the set decimal

    :param n: 1-D array vector
    :param decimal: Default 3
    :return: Angle between 2 vectors. Value in [0, 180]
    """
    n = np.around(n, decimal)
    return n


def get_angle(v1: FloatArray, v2: FloatArray) -> float:
    """Compute the angle between 2 vectors in degrees

    :param v1: 1-D array vector
    :param v2: 1-D array vector
    :return: Angle between 2 vectors. Value in [0, 180]
    """
    v1_norm, v2_norm = vector_norm(v1), vector_norm(v2)
    if not all([v1_norm, v2_norm]):
        return np.NAN

    cos = dot_product(v1, v2) / (v1_norm * v2_norm)
    return np.rad2deg(np.arccos(cos))


def get_dihedral_angle(v1: FloatArray, v2: FloatArray, v3: FloatArray) -> float:
    """Compute dihedral angle between 2 sub-planes formed by 3 vectors. Return degrees belongs to (-180, 180].

    For points that are sequentially numbered and located at positions p1, p2, p3, p4, bond vectors are defined by
    v1=p2−p1, p2=p3−p2, and v3=p4−p3, as the input of this function.

    :param v1: 1-D array vector
    :param v2: 1-D array vector
    :param v3: 1-D array vector
    :return: Dihedral angle in degrees belongs to (-180, 180]
    """
    if not all([vector_norm(v) for v in (v1, v2, v3)]):
        return np.NAN

    norm_p1, norm_p2 = cross_product(v1, v2), cross_product(v2, v3)
    p1_x_p2 = cross_product(norm_p1, norm_p2)

    y = dot_product(p1_x_p2, v2) * (1.0 / vector_norm(v2))
    x = dot_product(norm_p1, norm_p2)

    return np.degrees(np.arctan2(y, x))


def hamming_distance(str1: str, str2: str) -> Tuple[int, Tuple[str, str]]:
    """Compute hamming distance between 2 strings

    Raises :class:`ValueError` if 2 input have different length.

    :param str1: First input string
    :param str2: Second input string
    :return: Number of different characters in the strings, with 2 strings consisting them
    """
    count = 0
    length = len(str1)
    chain_dif1 = []
    chain_dif2 = []
    if length != len(str2):
        raise ValueError("Input lengths not equal")

    for i in range(length):
        if str1[i] != str2[i]:
            count += 1
            chain_dif1.append(str1[i])
            chain_dif2.append(str2[i])
    return count, (chain_dif1, chain_dif2)


def random_matrix(shape: Union[list, tuple], epsilon: float) -> FloatArray:
    """

    :param shape: First input string
    :param epsilon: Second input string
    :return:
    """
    u = np.random.normal(0, 1, shape)
    norm = np.sum(u**2, axis=1) ** 0.5
    r = np.random.uniform(0, epsilon, shape[0])
    x = (
        r.reshape(-1, 1) * u / norm.reshape(-1, 1)
    )  # scale vectors to restrict their distance within epsilon
    return x


def get_third_side_length(a: float, b: float, angle_deg: float) -> float:
    """Law of cosines

    :param a: side 1
    :param b: side 2
    :param angle_deg: side angle
    :return: the third side
    """
    angle_rad = np.deg2rad(angle_deg)
    c_squared = a**2 + b**2 - 2 * a * b * np.cos(angle_rad)
    c_squared = max(c_squared, 0.0)
    return np.sqrt(c_squared)


def get_RMSD(x: FloatArray, y: FloatArray, n_dim: int = 3) -> float:
    """Computes RMSD between 2 point sets.

    :param x: Point set 1.
    :param y: Point set 2.
    :param n_dim: Dimension of points. Default by 3 for 3D structures.
    :return: RMSD value.
    """
    N = np.size(x) / n_dim  # number of points
    return np.sqrt(np.sum((x - y) ** 2) / N)
