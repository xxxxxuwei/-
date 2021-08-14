# This file is part of pyTSEB for estimating the resistances to momentum and heat transport 
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
This module includes functions for calculating the resistances for
heat and momentum trasnport for both One- and Two-Source Energy Balance models.
Additional functions needed in are imported from the following packages
该模块包括计算一源和双源能量平衡模型的热和动量传输阻力的功能。
* :doc:`meteoUtils` for the estimation of meteorological variables.
* :doc:`MOsimilarity` for the estimation of the Monin-Obukhov length and stability functions.

PACKAGE CONTENTS
================
Resistances阻抗
-----------
* :func:`CalcR_A` Aerodynamic resistance.空气动力阻抗
* :func:`CalcR_S_Kustas` [Kustas1999]_ soil resistance.土壤阻抗
* :func:`CalcR_X_Norman` [Norman1995]_ canopy boundary layer resistance.冠层边界层阻抗

Stomatal conductance气孔导度
--------------------
* :func:`CalcStomatalConductanceTSEB` TSEB stomatal conductance.TSEB气孔导度
* :func:`CalcStomatalConductanceOSEB` OSEB stomatal conductance.OSEB气孔导度
* :func:`CalcCoef_m2mmol` Conversion factor from stomatal conductance from m s-1 to mmol m-2 s-1.气孔导度的转换因子

Estimation of roughness估计粗糙度
-----------------------
* :func:`CalcD_0` Zero-plane displacement height.零平面位移高度
* :func:`CalcRoughness` Roughness for different land cover types.不同地表覆盖类型的粗糙度
* :func:`CalcZ_0M` Aerodynamic roughness lenght.空气动力学粗糙度长度
* :func:`Raupach` Roughness and displacement height factors for discontinuous canopies.不连续冠层的粗糙度和位移高度因子

"""

#==============================================================================
# List of constants used in TSEB model and sub-routines   
#==============================================================================
# Land Cover Classes
CROP=11
GRASS=2
SHRUB=5
CONIFER=4
BROADLEAVED=3
# von Karman's constant
k=0.4 
#acceleration of gravity (m s-2)
gravity=9.8
# Universal gas constant (kPa m3 mol-1 K-1)
R_u=0.0083144

KN_b = 0.012 # Value propoesd in Kustas et al 1999
KN_c = 0.0025 # Coefficient from Norman et al. 1995
KN_C_dash = 90.0 # value proposed in Norman et al. 1995

import src.MOsimilarity as MO
import src.meteoUtils as met

def CalcD_0 (h_C):
    ''' Zero-plane displacement height零平面位移高度

    Calculates the zero-plane displacement height based on a 
    fixed ratio of canopy height.
    基于冠层高度的固定比率计算零平面位移高度
    Parameters
    ----------
    h_C : float
        canopy height (m).冠层高度
    
    Returns
    -------
    d_0 : float
        zero-plane displacement height (m).零平面位移高度'''
    
    d_0 = h_C * 0.65
    return d_0

def CalcRoughness (LAI, hc,wc=1,landcover=11):
    ''' Surface roughness and zero displacement height for different vegetated surfaces.
    不同植被表面粗糙度和零位移高度
    Calculates the roughness using different approaches depending we are dealing with
    crops or grasses (fixed ratio of canopy height) or shrubs and forests,depending of LAI
    and canopy shape, after [Schaudt2000]_
    根据我们处理的作物或草(冠层高度的固定比率)或灌木和森林，使用不同的方法计算粗糙度，这取决于叶面积指数和冠层形状
    Parameters
    ----------
    LAI : float
        Leaf (Plant) Area Index.叶面积指数
    hc : float
        Canopy height (m)冠层高度
    wc : float, optional
        Canopy height to width ratio.冠层高宽比
    landcover : int, optional
        landcover type, use 11 for crops, 2 for grass, 5 for shrubs,
        4 for conifer forests and 3 for broadleaved forests.
    
    Returns
    -------
    z_0M : float
        aerodynamic roughness length for momentum trasport (m).用于动量传输的空气动力粗糙度长度
    d : float
        Zero-plane displacement height (m).零平面位移高度
    
    References
    ----------
    .. [Schaudt2000] K.J Schaudt, R.E Dickinson, An approach to deriving roughness length
        and zero-plane displacement height from satellite data, prototyped with BOREAS data,
        Agricultural and Forest Meteorology, Volume 104, Issue 2, 8 August 2000, Pages 143-155,
        http://dx.doi.org/10.1016/S0168-1923(00)00153-2.
    '''

    from math import exp,pi
    #Needleleaf canopies
    if landcover == CONIFER:
        fc=1.-exp(-0.5*LAI)
        lambda_=(2./pi)*fc*wc
        #Calculation of the Raupach (1994) formulae
        z0M_factor,d_factor=Raupach(lambda_)
    #Broadleaved canopies
    elif landcover == BROADLEAVED:
        fc=1.-exp(-LAI)
        lambda_=fc*wc
        z0M_factor,d_factor=Raupach(lambda_)
    #Shrublands
    elif landcover == SHRUB:
        fc=1.-exp(-0.5*LAI)
        lambda_=fc*wc
        z0M_factor,d_factor=Raupach(lambda_)
    else:
        z0M_factor=0.125
        d_factor=0.65

    #Calculation of correction factors from  Lindroth
    if LAI <= 0:
        fz= 1.0
        fd= 1.0
    elif LAI < 0.8775:
        fz= 0.3299*LAI**1.5+2.1713
        fd=1.-0.3991*exp(-0.1779*LAI)
    else:
        fz=1.6771*exp(-0.1717*LAI)+1.
        fd=1.-0.3991*exp(-0.1779*LAI)
    #Application of the correction factors to roughness and displacement height
    z0M_factor=z0M_factor*fz
    d_factor=d_factor*fd
    if landcover == CROP or landcover == GRASS:
        z0M_factor=1./8.
        d_factor=0.65
    #Calculation of rouhgness length
    z_0M=z0M_factor*hc
    #Calculation of zero plane displacement height
    d=d_factor*hc
    return z_0M, d

def CalcR_A (z_T, ustar, L, d_0, z_0H, useRi=False, z_star=False):
    ''' Estimates the aerodynamic resistance to heat transport based on the
    MO similarity theory.
    基于MO相似理论估算了热输运的气动阻力。
    Parameters
    ----------
    z_T : float
        air temperature measurement height (m).气温测量高度
    ustar : float
        friction velocity (m s-1).摩擦速度
    L : float
        Monin Obukhov Length or Richardson number for stability (see useRi variable).稳定性的长度或理查森数
    d_0 : float
        zero-plane displacement height (m).零平面位移高度
    z_0M : float
        aerodynamic roughness length for momentum trasport (m).用于动量传输的空气动力粗糙度长度
    z_0H : float
        aerodynamic roughness length for heat trasport (m).用于热传输的空气动力粗糙度长度
    useRi : bool, optional
        boolean variable to use Richardsond number instead of the MO length.布尔变量，使用richardsonnumber代替MO的长度
    z_star : float or None, optional
        height of the roughness sublayer (RSL)粗糙度亚层高度(RSL), optional, if used and zstar>0 
        the adiabatic correction in the RSL will be computed.
    
    Returns
    -------
    R_A : float
        aerodyamic resistance to heat transport in the surface layer (s m-1).
        在表面层中对热传输的空气动力阻力
    References
    ----------        
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293, http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    '''

    from math import log
    if ustar==0:return float('inf')
    Psi_H_star=0.0
    R_A_log = log( (z_T - d_0) / z_0H)
    if (useRi):
        # use the approximation Ri ~ (z-d_0)./L from end of section 2.2 from Norman et. al., 2000 (DTD paper) 				
        Psi_H = MO.CalcPsi_H(L)
        Psi_H0 = MO.CalcPsi_H(L/(z_T - d_0)*z_0H)
    else:
        #if L -> infinity, z./L-> 0 and there is neutral atmospheric stability 
        Psi_H = 0.0
        Psi_H0 = 0.0
        #other atmospheric conditions
        if L == 0.0:L=1e-36
        Psi_H = MO.CalcPsi_H((z_T - d_0)/L)
        Psi_H0 = MO.CalcPsi_H(z_0H/L)
    if z_star>0 and z_T<=z_star:
        Psi_H_star=MO.CalcPsi_H_star(z_T, L, d_0,z_0H,z_star)
    R_A =  (R_A_log - Psi_H + Psi_H0+Psi_H_star) /(ustar * k)
    return R_A
   
def CalcR_S_Kustas (u_S, deltaT):
    ''' Aerodynamic resistance at the  soil boundary layer.
    土壤边界层空气动力阻力
    Estimates the aerodynamic resistance at the  soil boundary layer based on the
    original equations in TSEB [Kustas1999]_.根据TSEB的原始方程估计了土壤边界层的气动阻力

    Parameters
    ----------
    u_S : float
        wind speed at the soil boundary layer (m s-1).土壤边界层风速
    deltaT : float
        Surface to air temperature gradient (K).地面到空气的温度梯度
    
    Returns
    -------
    R_S : float
        Aerodynamic resistance at the  soil boundary layer (s m-1).
        土壤边界层空气动力阻力
    References
    ----------
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''
    if u_S==0:return float('inf')

    if deltaT<0.0:
        deltaT = 0.0
    R_S = 1.0/ (KN_c* deltaT**(1.0/3.0)+ KN_b * u_S)
    return R_S

def CalcR_X_Norman(LAI, leaf_width, u_d_zm):
    ''' Estimates aerodynamic resistance at the canopy boundary layer.
        估计冠层边界层的空气动力阻力
    Estimates the aerodynamic resistance at the  soil boundary layer based on the
    original equations in TSEB [Norman1995]_.根据TSEB的原始方程估计了冠层边界层的空气动力阻力

    Parameters
    ----------
    F : float
        local Leaf Area Index.局部叶面积指数
    leaf_width : float
        efective leaf width size (m).有效叶宽大小
    u_d_zm : float
        wind speed at the height of momomentum source-sink. .动量源-汇聚高度的风速
    
    Returns
    -------
    R_x : float
        Aerodynamic resistance at the canopy boundary layer (s m-1).
        冠层边界层的空气动力阻力
    References
    ----------    
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293, http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    '''

    if u_d_zm==0:return float('inf')
    #C_dash = 130.0 # Original value proposed by McNaughton & Van der Hurk 1995
    C_dash_F = KN_C_dash/LAI
    R_x = C_dash_F*(leaf_width/u_d_zm)**0.5
    return R_x
  
def CalcZ_0H (z_0M,kB=2):
    '''Estimate the aerodynamic routhness length for heat trasport.
        估计热传输的空气动力粗糙度长度
    Parameters
    ----------
    z_0M : float
        aerodynamic roughness length for momentum transport (m).用于动量传输的空气动力粗糙度长度
    kB : float
        kB parameter, default = 0.
    
    Returns
    -------
    z_0H : float
        aerodynamic roughness length for momentum transport (m).用于热传输的空气动力粗糙度长度
    
    References
    ----------
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293, http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    '''

    from math import exp
    z_OH = z_0M/exp(kB)
    return z_OH
    
def CalcZ_0M (h_C):
    ''' Aerodynamic roughness lenght.空气动力学粗糙度长度

    Estimates the aerodynamic roughness length for momentum trasport 
    as a ratio of canopy height.
    以冠层高度的比率估计动量传输的空气动力学粗糙度长度
    Parameters
    ----------
    h_C : float
        Canopy height (m).冠层高度
    
    Returns
    -------
    z_0M : float
        aerodynamic roughness length for momentum transport (m).用于动量传输的空气动力粗糙度长度'''

    z_OM = h_C * 0.125
    return z_OM

def Raupach(lambda_):
    '''Roughness and displacement height factors for discontinuous canopies
        不连续冠层的粗糙度和位移高度因子
    Estimated based on the frontal canopy leaf area, based on Raupack 1994 model,基于Raupack 1994模型，基于额冠层叶面积估算
    after [Schaudt2000]_
    
    Parameters
    ----------
    lambda_ : float
        roughness desnsity or frontal area index.
        粗糙度密度或锋面面积指数
    Returns
    -------
    z0M_factor : float
        height ratio of roughness length for momentum transport动量传输的粗糙度长度的高度比
    d_factor : float
        height ratio of zero-plane displacement height零平面位移高度的高度比

    References
    ----------
    .. [Schaudt2000] K.J Schaudt, R.E Dickinson, An approach to deriving roughness length
        and zero-plane displacement height from satellite data, prototyped with BOREAS data,
        Agricultural and Forest Meteorology, Volume 104, Issue 2, 8 August 2000, Pages 143-155,
        http://dx.doi.org/10.1016/S0168-1923(00)00153-2.

    '''   

    from math import exp,sqrt
    z0M_factor=0.125
    d_factor=0.65
    # Calculation of the Raupach (1994) formulae
    if lambda_ > 0.152:
        z0M_factor=(0.0537/(lambda_**0.510))*(1.-exp(-10.9*lambda_**0.874))+0.00368
    else:
        z0M_factor=5.86*exp(-10.9*lambda_**1.12)*lambda_**1.33+0.000860
    if lambda_ > 0: d_factor=1.-(1.-exp(-sqrt(15.0*lambda_)))/sqrt(15.0*lambda_)
    return z0M_factor,d_factor
