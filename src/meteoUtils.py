# This file is part of pyTSEB for calculating the meteorolgical variables
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

'''
Created on Apr 6 2015
@author: Hector Nieto (hnieto@ias.csic.es)

Modified on feb 3 2016
@author: Hector Nieto (hnieto@ias.csic.es)

DESCRIPTION
===========
This package contains functions for estimating meteorological variables needed
in resistance energy balance models.
该软件包包含用于估计阻力能量平衡模型所需的气象变量的功能

PACKAGE CONTENTS
================
* :func:`CalcC_p` Heat capacity of air at constant pressure.空气在恒定压力下的热容
* :func:`CalcLambda(Ta_K)` Latent heat of vaporization.汽化潜热
* :func:`CalcPressure` Barometric pressure.气压
* :func:`CalcPsicr` Psicrometric constant.Psicrometric常数
* :func:`CalcRho` Density of air.空气密度
* :func:`CalcStephanBoltzmann` Stephan-Boltzmann law for blackbody radiation emission.黑体辐射辐射的Stephan-Boltzmann定律
* :func:`CalcTheta_s` Sun Zenith Angle.太阳天顶角
* :func:`Get_SunAngles` Sun Zenith and Azimuth Angles.太阳天顶和方位角
* :func:`CalcVaporPressure` Saturation water vapour pressure.饱和水汽压
* :func:`CalcDeltaVaporPressure` Slope of saturation water vapour pressure.饱和水汽压的斜率
'''

#==============================================================================
# List of constants used in Meteorological computations  
#==============================================================================
#Stephan Boltzmann constant (W m-2 K-4)斯蒂芬玻尔兹曼常数
sb=5.670373e-8 
# heat capacity of dry air at constant pressure (J kg-1 K-1)干燥空气在恒定压力下的热容量
c_pd=1003.5 
# heat capacity of water vapour at constant pressure (J kg-1 K-1)恒压下水蒸气的热容
c_pv=1865
#ratio of the molecular weight of water vapor to dry air水蒸气与干燥空气的相对分子质量之比
epsilon=0.622 
#Psicrometric Constant kPa K-1   Psicrometric常数
psicr=0.0658    
# gas constant for dry air, J/(kg*degK) 干燥空气的气体常数
R_d=287.04


def CalcC_p (p, ea):
    ''' Calculates the heat capacity of air at constant pressure.
    计算空气在定压下的热容
    Parameters
    ----------
    p : float
        total air pressure (dry air + water vapour) (mb).总气压(干空气+水蒸气)
    ea : float
        water vapor pressure at reference height above canopy (mb).冠层上方参考高度的水汽压
    
    Returns
    -------
    c_p : heat capacity of (moist) air at constant pressure (J kg-1 K-1).(湿)空气在恒定压力下的热容量
    
    References
    ----------
    based on equation (6.1) from Maarten Ambaum (2010): 
    Thermal Physics of the Atmosphere (pp 109).'''

    # first calculate specific humidity, rearanged eq (5.22) from Maarten Ambaum (2010), (pp 100)
    q = epsilon * ea / (p + (epsilon - 1.0)*ea)
    # then the heat capacity of (moist) air
    c_p = (1.0-q)*c_pd + q*c_pv
    return c_p
    
def CalcLambda(Ta_K):
    '''Calculates the latent heat of vaporization.
    计算蒸发潜热
    Parameters
    ----------
    Ta_K : float
        Air temperature (Kelvin).空气温度
    
    Returns
    -------
    Lambda : float
        Latent heat of vaporisation (MJ kg-1).蒸发潜热
    
    References
    ----------
    based on Eq. 3-1 Allen FAO98 '''
    Lambda = 2.501 - (2.361e-3* (Ta_K-273.15) ) 
    return Lambda

def CalcPressure(z):
    ''' Calculates the barometric pressure above sea level.
    计算海平面以上的气压。
    Parameters
    ----------
    z: float
        height above sea level (m).海拔高度
    
    Returns
    -------
    p: float
        air pressure (mb).气压'''
    
    p=1013.25*(1.0-2.225577e-5*z)**5.25588
    return p

def CalcPsicr(p,Lambda):
    ''' Calculates the psicrometric constant.计算psicrometric常数
     
    Parameters
    ----------
    p : float
        atmopheric pressure (mb).气压
    Lambda : float
        latent heat of vaporzation (MJ kg-1).蒸发潜热
    
    Returns
    -------
    psicr : float
        Psicrometric constant (mb C-1).psicrometric常数'''
    
    psicr=0.00163*p/(Lambda)
    return psicr
    
def CalcRho (p, ea, Ta_K):
    '''Calculates the density of air.计算空气的密度。
    
    Parameters
    ----------
    p : float
        total air pressure (dry air + water vapour) (mb).总气压(干空气+水蒸气)
    ea : float
        water vapor pressure at reference height above canopy (mb).冠层上方参考高度的水汽压
    Ta_K : float
        air temperature at reference height (Kelvin).参考高度的气温
    
    Returns
    -------
    rho : float
        density of air (kg m-3).空气的密度

    References
    ----------
    based on equation (2.6) from Brutsaert (2005): Hydrology - An Introduction (pp 25).'''
    
    # p is multiplied by 100 to convert from mb to Pascals
    rho = ((p*100.0)/(R_d * Ta_K )) * (1.0 - (1.0 - epsilon) * ea / p)
    return rho

def CalcStephanBoltzmann(T_K):
    '''Calculates the total energy radiated by a blackbody.计算黑体辐射的总能量
    
    Parameters
    ----------
    T_K : float
        body temperature (Kelvin)体温
    
    Returns
    -------
    M : float
        Emitted radiance (W m-2)发射亮度'''
    
    M=sb*T_K**4
    return M
    
def CalcTheta_s (xlat,xlong,stdlng,doy,year,ftime):
    '''Calculates the Sun Zenith Angle (SZA).
    计算太阳天顶角
    Parameters
    ----------
    xlat : float
        latitude of the site (degrees).网站纬度
    xlong : float
        longitude of the site (degrees).网站的经度
    stdlng : float
        central longitude of the time zone of the site (degrees).站点时区的中央经度
    doy : float
        day of year of measurement (1-366).测量年的一天
    year : float
        year of measurement .年的测量
    ftime : float
        time of measurement (decimal hours).时间的测量
    
    Returns
    -------
    theta_s : float
        Sun Zenith Angle (degrees).太阳天顶角
    
    References
    ----------
    Adpopted from Martha Anderson's fortran code for ALEXI which in turn was based on Cupid.'''
    
    from math import pi, radians,sin, cos, asin, acos, degrees
    pid180 = pi/180
    pid2 = pi/2.0
    # Latitude computations
    xlat=radians(xlat)
    sinlat=sin(xlat)
    coslat=cos(xlat)
    # Declination computations
    kday=(year-1977.0)*365.0+doy+28123.0
    xm=radians(-1.0+0.9856*kday)
    delnu=2.0*0.01674*sin(xm)+1.25*0.01674*0.01674*sin(2.0*xm)
    slong=radians((-79.8280+0.9856479*kday))+delnu
    decmax=sin(radians(23.44))
    decl=asin(decmax*sin(slong))
    sindec=sin(decl)
    cosdec=cos(decl)
    eqtm=9.4564*sin(2.0*slong)/cosdec-4.0*delnu/pid180
    eqtm=eqtm/60.0
    # Get sun zenith angle
    timsun=ftime    #MODIS time is already solar time
    hrang=(timsun-12.0)*pid2/6.0
    theta_s = acos(sinlat*sindec+coslat*cosdec*cos(hrang))
    # if the sun is below the horizon just set it slightly above horizon
    theta_s = min(theta_s, pid2-0.0000001)
    theta_s=degrees(theta_s)
    return theta_s

def Get_SunAngles(lat,lon,stdlon,doy,ftime):
    '''Calculates the Sun Zenith and Azimuth Angles (SZA & SAA).
    计算太阳天顶和方位角
    Parameters
    ----------
    lat : float
        latitude of the site (degrees).网站纬度
    long : float
        longitude of the site (degrees).网站的经度
    stdlng : float
        central longitude of the time zone of the site (degrees).站点时区的中央经度
    doy : float
        day of year of measurement (1-366).测量年的一天
    ftime : float
        time of measurement (decimal hours).时间的测量
    
    Returns
    -------
    sza : float
        Sun Zenith Angle (degrees).太阳天顶角
    saa : float
        Sun Azimuth Angle (degrees).太阳方位角
     
    '''
  
    from math import pi, sin, cos, asin, acos, radians, degrees
    # Calculate declination
    declination=0.409*sin((2.0*pi*doy/365.0)-1.39)
    EOT=0.258*cos(declination)-7.416*sin(declination)-3.648*cos(2.0*declination)-9.228*sin(2.0*declination)
    LC=(stdlon-lon)/15.
    time_corr=(-EOT/60.)+LC
    solar_time=ftime-time_corr
    # Get the hour angle
    w=(solar_time-12.0)*15.
    # Get solar elevation angle
    sin_thetha=cos(radians(w))*cos(declination)*cos(radians(lat))+sin(declination)*sin(radians(lat))
    sun_elev=asin(sin_thetha)
    # Get solar zenith angle
    sza=pi/2.0-sun_elev
    sza=degrees(sza)
    # Get solar azimuth angle
    cos_phi=(sin(declination)*cos(radians(lat))-cos(radians(w))*cos(declination)*sin(radians(lat)))/cos(sun_elev)
    if w <= 0.0:
        saa=degrees(acos(cos_phi))
    else:
        saa=360.-degrees(acos(cos_phi))
    return sza,saa

def CalcVaporPressure(T_K):
    '''Calculate the saturation water vapour pressure.计算饱和水汽压
    
    Parameters
    ----------
    T_K : float
        temperature (K).温度
    
    Returns
    -------
    ea : float
        saturation water vapour pressure (mb).饱和水汽压'''
   
    from math import exp
    T_C=T_K-273.15
    ea= 6.112 * exp((17.67*T_C)/(T_C + 243.5))
    return ea
    
def CalcDeltaVaporPressure(T_K):
    '''Calculate the slope of saturation water vapour pressure.
    计算饱和水汽压的斜率
    Parameters
    ----------
    T_K : float
        temperature (K).温度
    
    Returns
    -------
    s : float
        slope of the saturation water vapour pressure (kPa K-1)饱和水汽压的斜率'''
   
    from math import exp
    T_C=T_K-273.15
    s= 4098.0 * (0.6108*exp(17.27*T_C/(T_C+237.3)))/((T_C+237.3)**2)
    return s
