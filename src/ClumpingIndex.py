# This file is part of pyTSEB for calculating the canopy clumping index
# Copyright 2016 Hector Nieto and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Created on Apr 6 2015
@author: Hector Nieto (hnieto@ias.csic.es)

Modified on Jan 27 2016
@author: Hector Nieto (hnieto@ias.csic.es)

DESCRIPTION
===========
Routines for calculating the clumping index for both randomly placed canopies and
structured row crops such as vineyards.
计算随机放置的树冠和结构行作物(如葡萄园)的聚集指数的例行程序。
PACKAGE CONTENTS
================
* :func:`CalcOmega0_Kustas` Nadir viewing clmping factor.最低点观测钳位因子
* :func:`CalcOmega_Kustas` Clumping index at an incidence angle.在入射角处的聚集指数
"""
def CalcOmega0_Kustas(LAI, f_C,x_LAD=1,isLAIeff=True):
    ''' Nadir viewing clmping factor最低点观测钳位因子

    Estimates the clumping factor forcing equal gap fraction between the real canopy
    and the homogeneous case, after [Kustas1999]_.
    估计了迫使真实冠层与均匀冠层之间的间隙分数相等的聚集因子     
    Parameters
    ----------  
    LAI : float
        Leaf Area Index, it can be either the effective LAI or the real LAI  叶面积指数，可以是有效叶面积指数，也可以是真实叶面积指数
        , default input LAI is effective. 默认输入有效叶面积指数
    f_C : float
        Apparent fractional cover, estimated from large gaps, means that
        are still gaps within the canopy to be quantified.明显的部分覆盖，从大间隙估计，意味着仍是树冠内的间隙需要量化。
    x_LAD : float, optional
        Chi parameter for the ellipsoildal Leaf Angle Distribution function of 参数为椭球叶角分布函数
        [Campbell1988]_ [default=1, spherical LIDF].
    isLAIeff :  bool, optional
        Defines whether the input LAI is effective or local.
    定义输入叶面积指数是有效的还是局部的
    Returns
    -------
    omega0 : float
        clumping index at nadir.在最低点的聚集指数

    References
    ----------
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
 '''
    
    from math import log,exp, sqrt,radians, tan
    
    theta=0.0
    theta=radians(theta)    
    # Estimate the beam extinction coefficient based on a ellipsoidal LAD function
    # Eq. 15.4 of Campbell and Norman (1998)
    K_be=sqrt(x_LAD**2+tan(theta)**2)/(x_LAD+1.774*(x_LAD+1.182)**-0.733)
    if isLAIeff:
        F=LAI/f_C
    else: # The input LAI is actually the real LAI
        F=float(LAI)
    # Calculate the gap fraction of our canopy
    trans = f_C*exp(-K_be * F)+(1.0-f_C)
    if trans<=0:
        trans=1e-36
    # and then the nadir clumping factor
    omega0 = -log(trans)/(LAI*K_be)
    return omega0

def CalcOmega_Kustas(omega0,theta,wc=1):
    ''' Clumping index at an incidence angle.在入射角处的聚集指数

    Estimates the clumping index for a given incidence angle assuming randomnly placed canopies.
    估计一个给定入射角假设随机放置的冠层的聚集指数
    Parameters
    ----------
    omega0 : float
        clumping index at nadir在最低点的聚集指数, estimated for instance by :func:`CalcOmega0_Kustas`.
    theta : float
        incidence angle (degrees).入射角
    wc :  float, optional
        canopy witdth to height ratio冠层高度比, [default = 1].

    Returns
    -------
    Omega : float
        Clumping index at an incidenc angle.在入射角的聚集指数

    References
    ----------
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''
    
    from math import exp,radians
    wc=1.0/wc
    omega = omega0 / (omega0 + (1.0 - omega0) * exp(-2.2 * (radians(theta))**(3.8 - 0.46 * wc)))
    return omega
