# This file is part of pyTSEB for processes related to the Monin-Obukhov Similarity Theory 这个文件是pyTSEB的一部分，为了与莫宁-奥布霍夫长度相似理论相关的进程
# Copyright 2016 Hector Nieto and contributors listed in the README.md file. 版权所有2016 Hector Nieto 和 贡献者列出在README.md文件
# This program is free software: you can redistribute it and/or modify 这个程序是自由软件:你可以重新区分和/或修改它
# it under the terms of the GNU Lesser General Public License as published by 它在GNU Lesser通用公共许可证的条款下被发布
# the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
#自由软件基金会，或者授权的第三版，或者(由您选择)任何以后的版本。
# This program is distributed in the hope that it will be useful,发布这个程序是希望它能有用，
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#但没有任何保证;甚至没有暗示性的保证，针对特定目的的适销性或适用性。依据 GNU Lesser General Public License获取更多细节。
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#您应该已经收到GNU Lesser通用公共许可证的副本还有这个程序。如果不是，请参见<http://www.gnu.org/licenses/>。
"""
Created on Apr 6 2015
@author: Hector Nieto (hnieto@ias.csic.es)

Modified on Jan 27 2016 #修改
@author: Hector Nieto (hnieto@ias.csic.es)

DESCRIPTION 这个库的用途
===========
This package contains the main routines for estimating variables related to the Monin-Obukhov (MO) Similarity Theory,
such as  MO length, adiabatic correctors for heat and momentum transport and wind profile through a canopy. It
requires the following package.
这个库包含了估计与Monin-Obukhov (MO莫宁-奥布霍夫长度)相似理论相关的变量的主要程序，
例如MO长度、热和动量传输的绝热校正器以及通过冠层的风廓线。
它需要以下软件包。
* :doc:`meteoUtils` for the estimation of meteorological variables.
'meteoUtils'用于估计气象变量。

PACKAGE CONTENTS 这个库可以算哪些参数
================
* :func:`CalcL` Monin-Obukhov length. Monin-Obukhov长度。
* :func:`CalcRichardson` Richardson number. 理查森数:在大气上，理查森数表示大气静力稳定度与垂直风切变的比值。
* :func:`CalcU_star` Friction velocity.摩擦速度

Stability correction functions稳定性校正函数
------------------------------
* :func:`CalcPsi_H` Adiabatic correction factor for heat transport.热传输的绝热修正因子
* :func:`CalcPsi_M` Adiabatic correction factor for momentum transport.动量传输的绝热修正因子
* :func:`CalcPsi_M_B92` Adiabatic correction factor for momentum transport [Brutsaert1992]_.动量输运的绝热修正因子[Brutsaert1992]

Wind profile functions风速剖面函数/风廓线函数
----------------------
* :func:`CalcU_C` [Norman1995]_ canopy wind speed.冠层风速
* :func:`CalcU_Goudriaan` [Goudriaan1977]_ wind speed profile below the canopy.冠层下的风速分布图
* :func:`CalcA_Goudriaan` [Goudriaan1977]_ wind attenuation coefficient below the canopy.冠层以下的风衰减系数

"""

#==============================================================================
# List of constants used in MO similarity MO相似性中使用的常量列表
#==============================================================================
# von Karman's constant冯卡曼常数
k=0.4 
#acceleration of gravity (m s-2)重力加速度
gravity=9.8
# Drag coefficient (Goudriaan 1977)阻力系数
c_d=0.2
# Relative turbulence intensity (Goudriaan 1977)相对湍流强度
i_w=0.5

import meteoUtils as met

def CalcL (ustar, Ta_K, rho, c_p, H, LE):
    '''Calculates the Monin-Obukhov length.计算Monin-Obukhov长度

    Parameters参数
    ----------
    ustar : float
        friction velocity (m s-1).摩擦速度
    Ta_K : float
        air temperature (Kelvin).气温
    rho : float
        air density (kg m-3).空气密度
    c_p : float
        Heat capacity of air at constant pressure (J kg-1 K-1).空气在恒定压力下的热容
    H : float
        sensible heat flux (W m-2).感热通量
    LE : float
        latent heat flux (W m-2).潜热通量
    
    Returns
    -------
    L : float
        Obukhov stability length (m).大气稳定度参数，输入m计算
    
    References 疑问！！！
    ----------
    .. [Brutsaert2005] Brutsaert, W. (2005). Hydrology: an introduction (Vol. 61, No. 8).
        Cambridge: Cambridge University Press.'''

    # first convert latent heat into rate of surface evaporation (kg m-2 s-1)
    Lambda = met.CalcLambda(Ta_K)*1e6 #in J kg-1
    E = LE / Lambda
    # Virtual sensible heat flux
    Hv=H+(0.61*Ta_K*c_p*E)
    if Hv!=0:
        L_const = k*gravity/Ta_K
        L = -ustar**3 / ( L_const*(Hv/(rho*c_p) ))
    else:
        L = float('inf')
    return L
    
def CalcPsi_H (zoL):
    ''' Calculates the adiabatic correction factor for heat transport.计算了热输运的绝热修正系数
    
    Parameters
    ----------
    zoL : float
        stability coefficient (unitless).稳定性系数
    
    Returns
    -------
    Psi_H : float
        adiabatic corrector factor fof heat transport (unitless).绝热校正因子fof热传输
        
    References
    ----------
    .. [Brutsaert2005] Brutsaert, W. (2005). Hydrology: an introduction (Vol. 61, No. 8).
        Cambridge: Cambridge University Press.
    '''
    
    from math import log
    Psi_H = 0.0
    #for stable and netural (zoL = 0 -> Psi_H = 0) conditions
    if zoL>=0.0:
        a = 6.1
        b = 2.5
        Psi_H = -a * log( zoL + (1.0 + zoL**b)**(1./b))
    # for unstable conditions
    else:
        y = -zoL
        c = 0.33		
        d = 0.057
        n = 0.78
        Psi_H = ((1.0-d)/n) * log((c + y**n)/c)
    return Psi_H

def CalcPsi_M (zoL):
    ''' Adiabatic correction factor for momentum transport.动量传输的绝热修正因子
    
    Parameters
    ----------
    zoL : float
        stability coefficient (unitless).稳定系数
    
    Returns
    -------
    Psi_M : float
        adiabatic corrector factor fof momentum transport (unitless).绝热修正因子fof(野外观测站)动量传输

    References
    ----------
    .. [Brutsaert2005] Brutsaert, W. (2005). Hydrology: an introduction (Vol. 61, No. 8).
        Cambridge: Cambridge University Press.
    '''

    from math import log, pi,atan
    Psi_M = 0.0
    # for stable and netural (zoL = 0 -> Psi_M = 0) conditions
    if zoL>=0.0:
        a = 6.1 
        b = 2.5
        Psi_M = -a * log( zoL + (1.0 + zoL**b)**(1.0/b))
    # for unstable conditions
    else:
        y = -zoL
        a = 0.33
        b = 0.41
        x = (y/a)**0.333333
        Psi_0 = -log(a) + 3**0.5*b*a**0.333333*pi/6.0
        y = min(y, b**-3)
        Psi_M = (log(a + y) - 3.0*b*y**0.333333 + (b*a**0.333333)/2.0 * log((1.0+x)**2/(1.0-x+x**2))+
            3.0**0.5*b*a**0.333333*atan((2.0*x-1.0)/3**0.5) + Psi_0)
    return Psi_M

def CalcPsi_M_B92 (zoL):
    ''' Adiabatic correction factor for momentum transport [Brutsaert1992]_动量传输的绝热修正因子
    
    Parameters
    ----------
    zoL : float
        stability coefficient (unitless).稳定系数
    
    Returns
    -------
    Psi_M : float
        adiabatic corrector factor fof momentum transport (unitless).绝热修正因子fof动量传输(无单位的)
    
    References参考文献
    ----------
    .. [Brutsaert1992] Brutsaert, W. (1992). Stability correction functions for the mean wind speed and 稳定校正函数的平均风速和不稳定表面层的温度
        temperature in the unstable surface layer. Geophysical research letters, 19(5), 469-472,地球物理研究快报
        http://dx.doi.org/10.1029/92GL00084.
    '''

    from math import log, pi
    Psi_M = 0.0
    # for stable and netural (zoL = 0 -> Psi_M = 0) conditions, Eq. 2.59 in Brutasert 2005  
    if zoL>=0.0:
        a = 6.1 
        b = 2.5
        Psi_M = -a * log( zoL + (1.0 + zoL**b)**(1.0/b))
    # for unstable conditions
    else:
        y=-zoL
        if y <0.0059:
            Psi_M=0.0
        elif y<=15.025:
            y_c=0.0059
            Psi_M = 1.47*log((0.28+y**0.75)/(0.28+(0.0059+y_c)**0.75))-1.29*(y**(1./3.)-(0.0059+y_c)**(1./3.))
        else:
            y_c=0.0059
            y=15.025
            Psi_M = 1.47*log((0.28+y**0.75)/(0.28+(0.0059+y_c)**0.75))-1.29*(y**(1./3.)-(0.0059+y_c)**(1./3.))

            
    return Psi_M

def CalcRichardson (u, z_u, d_0, T_R0, T_R1, T_A0, T_A1):
    '''Richardson number.理查森数

    Estimates the Bulk Richardson number for turbulence using time difference temperatures.
    利用时差温度估算湍流的批量理查森数。
    Parameters
    ----------
    u : float
        Wind speed (m s-1).风速
    z_u : float
        Wind speed measurement height (m).依据风速测量高度
    d_0 : float
        Zero-plane displacement height (m).零面位移高度
    T_R0 : float
        radiometric surface temperature at time 0 (K).时间0(k)时的辐射表面温度
    T_R1 : float
        radiometric surface temperature at time 1 (K).时间1(k)时的辐射表面温度
    T_A0 : float
        air temperature at time 0 (K).时间0(k)时空气温度(K)
    T_A1 : float
        air temperature at time 1 (K).时间1(K)的空气温度
    Returns
    -------
    Ri : float
        Richardson number.理查森数

    References
    ----------
    .. [Norman2000] Norman, J. M., W. P. Kustas, J. H. Prueger, and G. R. Diak (2000),
        Surface flux estimation using radiometric temperature: A dual-temperature-difference
        method to minimize measurement errors, Water Resour. Res., 36(8), 2263-2274,
        http://dx.doi.org/10.1029/2000WR900033.
    '''

    # See eq (2) from Louis 1979
    Ri = -(gravity * (z_u - d_0) / T_A1) * (((T_R1 - T_R0) - (T_A1 - T_A0)) / u**2)
    return Ri

def CalcU_C (u, h_C, d_0, z_0M, z_u,L):
    '''[Norman1995]_ wind speed at the canopy冠层处的风速
    
    Parameters
    ----------
    u : float
        Wind speed measured at heigth z_u (m s-1).高度z u (m s-1)测得的风速
    h_C : float
        canopy height (m).冠层高度
    d_0 : float
        zero-plane displacement height.零面位移高度
    z_0M : float
        aerodynamic roughness length for momentum transport (m).用于动量传输的空气动力学粗糙度
    z_u:  float
        Height of measurement of wind speeed.测量风速的高度
    L : float
     Monin Obukhov length.莫宁奥布霍夫长度
    
    Returns
    -------
    u_C : float
        wind speed at the canop interface (m s-1).冠层界面的风速
    
    References
    ----------
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293,
        http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    '''

    from math import log
    Psi_M= CalcPsi_M((h_C - d_0)/L)
    # calcualte u_C, wind speed at the top of (or above) the canopy
    u_C = u*log ((h_C  - d_0) / z_0M)/(log ((z_u  - d_0) / z_0M)- Psi_M)
    return u_C

def CalcU_Goudriaan (u_C, h_C, LAI, leaf_width, z):
    '''Estimates the wind speed at a given height below the canopy.
    估算树冠下给定高度下的风速
    Parameters
    ----------
    U_C : float
        Windspeed at the canopy interface (m s-1).冠层界面的风速
    h_C : float
        canopy height (m).树冠高度
    LAI : float
        Efective Leaf (Plant) Area Index.有效叶面积指数(植物)
    leaf_width : float
        effective leaf width size (m).有效叶宽大小(m)
    z : float
        heigh at which the windsped will be estimated.被估算风速处的高度
    
    Returns
    -------
    u_z : float
        wind speed at height z (m s-1).高度为z (m s-1)的风速
    
    References
    ----------
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293,
        http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    .. [Goudriaan1977] Goudriaan (1977) Crop micrometeorology: a simulation study
 '''

    from math import exp
    a=CalcA_Goudriaan (h_C,LAI,leaf_width) # extinction factor for wind speed
    u_z = u_C * exp(-a * (1.0 - (z/h_C))) # Eq. 4.48 in Goudriaan 1977
    return u_z

def CalcA_Goudriaan (h_C,LAI,leaf_width):
    ''' Estimates the extinction coefficient factor for wind speed估算消光系数因子的风速
    消光系数(exctinction coefficient)：是表征介质使电磁波衰减程度的物理量，它等于电磁波在介质中传播单位距离时,其强度由于吸收和散射作用而衰减的相对值。
    Parameters
    ----------    
    h_C : float
        canopy height (m)树冠高度
    LAI : float
        Efective Leaf (Plant) Area Index 有效叶面积指数(植物)
    leaf_width : float
        effective leaf width size (m)有效叶宽大小(m)

    Returns
    -------
    a : float
        exctinction coefficient for wind speed through the canopy 通过冠层的风速的消光系数
    
    References
    ----------
    .. [Goudriaan1977] Goudriaan (1977) Crop micrometeorology: a simulation study
    '''

    # Equation in Norman et al. 1995
    k3_prime=0.28
    a = k3_prime * LAI**(2./3.) * h_C**(1./3.) * leaf_width**(-1./3.)

    return a
    
def CalcU_star (u, z_u, L, d_0,z_0M, useRi=False):
    '''Friction velocity.摩擦速度
    Parameters
    ----------
    u : float
        wind speed above the surface (m s-1).地面以上风速(m s-1)
    z_u : float
        wind speed measurement height (m).风速测量高度
    L : float
        Obukhov stability length (m) or Richardson number, see useRi variable.奥布霍夫稳定长度(m)或理查森数，见useRi变量(见下面)
    d_0 : float
        zero-plane displacement height (m).零平面位移高度
    z_0M : float
        aerodynamic roughness length for momentum transport (m).用于动量传输的空气动力学粗糙度
    useRi : bool, optional 
        Use the Richardson number istead of MO length for adiabatic correction factors.使用理查森数代替MO长度的绝热修正因子
 
    References
    ----------
    .. [Brutsaert2005] Brutsaert, W. (2005). Hydrology: an introduction (Vol. 61, No. 8).
        Cambridge: Cambridge University Press.
    ''' 
    from math import log
    Psi_M=0.0
    Psi_M0=0.0
    if useRi: 
        #use the approximation Ri ~ (z-d_0)./L from end of section 2.	2 from Norman et. al., 2000 (DTD paper) 
        Psi_M = CalcPsi_M(L);
        Psi_M0 = CalcPsi_M(L/(z_u - d_0)*z_0M)
    else:
        #calculate correction factors in other conditions
        if L == 0.0: L=1e-36
        Psi_M= CalcPsi_M((z_u - d_0)/L)
        Psi_M0 = CalcPsi_M(z_0M/L)
    u_star = u * k / ( log( (z_u - d_0) / z_0M ) - Psi_M + Psi_M0)
    return u_star
