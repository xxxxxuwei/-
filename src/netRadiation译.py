# This file is part of pyTSEB for calculating the net radiation and its divergence
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

# 此文件是用于计算净辐射及其发散的 pyTSEB 的一部分
# 版权所有 2016 赫克托尼托和贡献者列在 README.md 文件。
#
# 此程序是免费软件：您可以重新分配和/或修改
#它根据GNU小一般公共许可证的条款，由
# 免费软件基础，许可证的第3版，或
#（在您的选项）任何较晚的版本。
#
# 此程序的分发希望它将是有用的，
#但没有任何保证：甚至没有隐含的保修
#特定目的的可商性或健身。 请参阅
# GNU 较小的一般公共许可证，了解更多详细信息。
#
#您应该已收到GNU小一般公共许可证的副本
#连同此程序。 如果没有，请参阅<http://www.gnu.org/licenses/>。
'''
Created on Apr 6 2015
@author: Hector Nieto (hnieto@ias.csic.es).

Modified on Jan 27 2016
@author: Hector Nieto (hnieto@ias.csic.es).

DESCRIPTION
===========
This package contains functions for estimating the net shortwave and longwave radiation
for soil and canopy layers. Additional packages needed are.

* :doc:`meteoUtils` for the estimation of meteorological variables.

PACKAGE CONTENTS
================
* :func:`CalcDifuseRatio` estimation of fraction of difuse shortwave radiation.
* :func:`CalcEmiss_atm` Atmospheric emissivity.
* :func:`CalcKbe_Campbell` Beam extinction coefficient.
* :func:`CalcLnKustas` Net longwave radiation for soil and canopy layers.
* :func:`CalcRnOSEB` Net radiation in a One Source Energy Balance model.
* :func:`CalcSnCampbell` Net shortwave radiation. 
'''

'''
此封装包含用于估计土壤和树冠层净短波和长波辐射的功能。
其他需要的封装是：
*:doc：用于估计气象变量的'气象图解'。

封装内容
================
*：func：'CalcDifuseRatio'估计短波辐射的分数。
*：func：'CalcEmiss_atm'大气的使者。
*：func：'CalcKbe_Campbell'消光系数。
*：func：'CalcLnKustas'土壤和树冠层的净长波辐射。
*：func：'CalcRnOSEB'单源能量平衡模型中的净辐射。
*：func：'CalcSnCampbell'净短波辐射。
'''

#==============================================================================
# List of constants used in the netRadiation Module净辐射模块中使用的常量列表
#==============================================================================
#Stephan Boltzmann constant (W m-2 K-4)
sb=5.670373e-8 

import meteoUtils as met

def CalcDifuseRatio(Sdn,sza,Wv=1,press=1013.25):
    '''Fraction of difuse shortwave radiation.漫射短波辐射的分数

    Partitions the incoming solar radiation into PAR and non-PR and
    diffuse and direct beam component of the solar spectrum.
    将入射的太阳辐射划分为PAR（光合有效辐射）和非PR，以及太阳光谱的散射和直接光束分量
    Parameters
    ----------
    Sdn : float
        Incoming shortwave radiation (W m-2).入射短波辐射
    sza : float
        Solar Zenith Angle (degrees).太阳天顶角
    Wv : float, optional
        Total column precipitable water vapour (g cm-2), default 1 g cm-2.
    press : float, optional
        atmospheric pressure (mb), default at sea level (1013mb).
        
    Returns
    -------
    difvis : float
        diffuse fraction in the visible region.在可见区域的扩散分数
    difnir : float
        diffuse fraction in the NIR region.在近红外区的扩散分数
    fvis : float
        fration of total visible radiation.总可见辐射的边缘
    fnir : float
        fraction of total NIR radiation.占近红外辐射总量的比例
    
    References
    ----------
    .. [Weiss1985] Weiss and Norman (1985) Partitioning solar radiation into direct and diffuse,
        visible and near-infrared components, Agricultural and Forest Meteorology,
        Volume 34, Issue 2, Pages 205-213,
        http://dx.doi.org/10.1016/0168-1923(85)90020-6.
    '''


    from math import radians, cos,exp, log10
    coszen=abs(cos(radians(sza)))
    #Calculate potential (clear-sky) visible and NIR solar components可见光和近红外太阳组件
    # Weiss & Norman 1985
    #Correct for curvature曲率纠正 of atmos in airmas (Kasten and Young,1989)
    if sza >90:
        difvis,difnir, fvis,fnir=[1.0,1.0, 0.4,0.6]
        return difvis,difnir, fvis,fnir
    else:
        #airmas=1.0/(abs(coszen)+0.50572*radians(96.07995-sza)**-1.6364)
        airmas=1.0/coszen
    #Visible PAR/NIR direct beam radiation可见光光合有效辐射/近红外直射辐射
    Rdirvis=600.*exp(-.185*airmas)*coszen                                   #Eq. 1
    Rdirvis=max(0,Rdirvis)
    w=1320.0*Wv*10**(-1.195+.4459*log10(airmas)-.0345*log10(airmas)**2)     #Eq. 6
    Rdirnir=(720.*exp(-0.06*(press/1313.25)*airmas)-w)*coszen                               #Eq. 4
    Rdirnir=max(0,Rdirnir)    
    # Potential diffuse radiation潜在漫射辐射
    Rdifvis=0.4*(600.0-Rdirvis)*coszen                                      #Eq. 3
    Rdifnir=0.6*(720.0-Rdirvis-w)*coszen                                    #Eq. 5
    Rdifvis=max(0,Rdifvis)
    Rdifnir=max(0,Rdifnir)      
    #Potential total solar radiation总太阳辐射
    potvis=Rdirvis+Rdifvis
    if potvis<=0:potvis=1e-6
    potnir=Rdirnir+Rdifnir
    if potnir<=0:potnir=1e-6
    fclear=Sdn/(potvis+potnir)
    fclear=min(1.0,fclear)
    #Partition SDN into VIS and NIR      
    fvis=potvis/(potvis+potnir)                                             #Eq. 7
    fnir=potnir/(potvis+potnir)                                            #Eq. 8
    fvis=max(0,fvis)
    fvis=min(1,fvis)
    fnir=1.0-fvis
    #Estimate direct beam and diffuse fractions in VIS and NIR wavebands
    ratiox=fclear
    if fclear > 0.9:ratiox=.9
    dirvis=(Rdirvis/potvis)*(1.-((.9-ratiox)/.7)**.6667)                    #Eq. 11
    if fclear > 0.88:ratiox=.88
    dirnir=(Rdirnir/potnir)*(1.-((.88-ratiox)/.68)**.6667)                  #Eq. 12
    dirvis=max(0.0,dirvis)
    dirnir=max(0.0,dirnir)
    dirvis=min(1.0,Rdirvis/potvis,dirvis)
    dirnir=min(1.0,Rdirnir/potnir,dirnir)
    if dirvis < 0.01 and dirnir > 0.01:dirvis=.011
    if dirnir < 0.01 and dirvis > 0.01:dirnir=.011
    difvis=1.0-dirvis
    difnir=1.0-dirnir
    return difvis,difnir, fvis,fnir

def CalcEmiss_atm(ea,Ta_K):
    '''Atmospheric emissivity大气发射率
    
    Estimates the effective atmospheric emissivity for clear sky.
    估计晴朗天空的大气有效发射率

    Parameters参数
    ----------
    ea : float
        atmospheric vapour pressure (mb).
        大气水汽压
    Ta_K : float
        air temperature (Kelvin).
        气温
    
    Returns
    -------
    emiss_air : float
        effective atmospheric emissivity.
        大气有效发射率

    References
    ----------    
    .. [Brutsaert1975] Brutsaert, W. (1975) On a derivable formula for long-wave radiation
        from clear skies, Water Resour. Res., 11(5), 742-744,
        htpp://dx.doi.org/10.1029/WR011i005p00742.'''

    emiss_air=1.24*(ea/Ta_K)**(1./7.)
    return emiss_air

def CalcKbe_Campbell(theta,x_LAD=1):
    ''' Beam extinction coefficient消光系数

    Calculates the beam extinction coefficient based on [Campbell1998]_ ellipsoidal
    leaf inclination distribution function.基于叶倾角分布椭球函数计算消光系数
    
    Parameters参数
    ----------
    theta : float
        incidence zenith angle (degrees).入射天顶角
    x_LAD : float, optional
        Chi parameter for the ellipsoidal Leaf Angle Distribution function, 叶倾角分布椭球函数的chi参数
        use x_LAD=1 for a spherical LAD.
    
    Returns
    -------
    K_be : float
        beam extinction coefficient.消光系数
    x_LAD: float, optional
        x parameter for the ellipsoidal Leaf Angle Distribution function, 叶倾角分布椭球函数的x参数
        use x_LAD=1 for a spherical LAD.
    
    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    '''
        
    from math import tan, radians,sqrt
    theta=radians(theta)
    K_be=sqrt(x_LAD**2+tan(theta)**2)/(x_LAD+1.774*(x_LAD+1.182)**-0.733)
    return K_be
    
def CalcLnKustas (T_C, T_S, Lsky, LAI, emisVeg, emisGrd,x_LAD=1):
    ''' Net longwave radiation for soil and canopy layers土壤和植被的长波净辐射

    Estimates the net longwave radiation for soil and canopy layers unisg based on equation 2a
    from [Kustas1999]_ and incorporated the effect of the Leaf Angle Distribution based on [Campbell1998]_
    
    Parameters
    ----------
    T_C : float
        Canopy temperature (K).
        冠层温度
    T_S : float
        Soil temperature (K).
        土壤温度
    Lsky : float
        Downwelling atmospheric longwave radiation (w m-2).
        下沉大气长波辐射
    LAI : float
        Effective Leaf (Plant) Area Index.
        有效叶（植物）面积指数
    emisVeg : float
        Broadband emissivity of vegetation cover.
        植被宽带发射率
    emisGrd : float
        Broadband emissivity of soil.
        土壤宽带发射率
    x_LAD: float, optional
        x parameter for the ellipsoidal Leaf Angle Distribution function, 
        use x_LAD=1 for a spherical LAD.

    Returns
    -------
    L_nC : float
        Net longwave radiation of canopy (W m-2).
        冠层净长波辐射
    L_nS : float
        Net longwave radiation of soil (W m-2).  
        土壤净长波辐射 

    References
    ----------
    .. [Kustas1999] Kustas and Norman (1999) Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''

    from math import exp,sqrt,log,cos,sin,radians
    # Integrate to get the diffuse transmitance积分以获得漫射透射率
    taud=0
    for angle in range(0,90,5):
        akd=CalcKbe_Campbell(angle,x_LAD)# Eq. 15.4
        taub=exp(-akd*LAI)
        taud = taud+taub*cos(radians(angle))*sin(radians(angle))*radians(5)
    taud=2.0*taud
    #D I F F U S E   C O M P O N E N T S
    #Diffuse light canopy reflection coefficients  for a deep canopy天篷漫射光冠层反射系数
    akd=-log(taud)/LAI
    ameanl=emisVeg    
    taudl=exp(-sqrt(ameanl)*akd*LAI)    #Eq 15.6
    # calculate long wave emissions from canopy, soil and sky计算冠层、土壤、大气长波辐射
    L_C = emisVeg*met.CalcStephanBoltzmann(T_C)
    L_S = emisGrd*met.CalcStephanBoltzmann(T_S)
    # calculate net longwave radiation divergence of the soil计算土壤净长波辐发射
    L_nS = taudl*Lsky + (1.0-taudl)*L_C - L_S
    L_nC = (1.0-taudl)*(Lsky + L_S - 2.0*L_C)
    return L_nC,L_nS

def CalcRnOSEB(Sdn,Lsky, T_R, emis, albedo):
    ''' Net radiation in a One Source Energy Balance model单源能量模型净辐射

    Estimates surface net radiation assuming a single `big leaf` layer.
    估算表面净辐射假设的一个‘大叶’层。

        
    Parameters参数
    ----------
    Sdn : float
        incoming shortwave radiation (W m-2).
        入射短波辐射
    Lsky : float
        Incoming longwave radiation (W m-2).   
        入射长波辐射
    T_R : float
        Radiometric surface temperature (K).
        辐射表面温度
    emis : float
        Broadband emissivity.
        宽带发射率
    albedoVeg : float
        Broadband short wave albedo.
        宽带短波反照率

    Returns
    -------
    R_s : float
        Net shortwave radiation (W m-2).
        净短波辐射
    R_l : float
        Net longwave radiation (W m-2).
        净长波辐射
    '''
    # outgoing shortwave radiation出射短波辐射
    R_sout = Sdn * albedo
    # outgoing long wave radiation出射长波辐射
    R_lout = emis * sb * (T_R)**4 + Lsky * (1.0 - emis)
    R_s = Sdn-R_sout
    R_l = Lsky-R_lout
    return R_s, R_l
    
def CalcSnCampbell (LAI, sza, Sdn_dir, Sdn_dif, fvis,fnir, rho_leaf_vis,
                    tau_leaf_vis,rho_leaf_nir, tau_leaf_nir, rsoilv, rsoiln,x_LAD=1):
    ''' Net shortwave radiation 净短波辐射

    Estimate net shorwave radiation for soil and canopy below a canopy using the [Campbell1998]_
    Radiative Transfer Model, and implemented in [Kustas1999]_
    估算土壤和冠层净短波辐射用树冠下辐射传输模型
    
    Parameters参数
    ----------
    LAI : float
        Effective Leaf (Plant) Area Index.
        有效叶（植物）面积指数
    sza : float
        Sun Zenith Angle (degrees).
        太阳高度角
    Sdn_dir : float
        Broadband incoming beam shortwave radiation (W m-2).
        宽带入射短波辐射
    Sdn_dir : float
        Broadband incoming beam shortwave radiation (W m-2).
    fvis : float
        fration of total visible radiation.
        总可见光辐射分数
    fnir : float
        fraction of total NIR radiation.
        近红外光辐射分数
    rho_leaf_vis : float
        Broadband leaf bihemispherical reflectance in the visible region (400-700nm).
        可见光域宽带叶双半球反射率
    tau_leaf_vis : float
        Broadband leaf bihemispherical transmittance in the visible region (400-700nm).
        可见光域宽带叶双半球透射比
    rho_leaf_nir : float
        Broadband leaf bihemispherical reflectance in the NIR region (700-2500nm).
        近红外域宽带叶双半球反射率
    tau_leaf_nir : float
        Broadband leaf bihemispherical transmittance in the NIR region (700-2500nm).
        近红外域宽带叶双半球透射比
    rsoilv : float
        Broadband soil bihemispherical reflectance in the visible region (400-700nm).
        可见光域宽带土壤双半球反射率
    rsoiln : float
        Broadband soil bihemispherical reflectance in the NIR region (700-2500nm).
        近红外域宽带土壤双半球反射率
    x_LAD : float,  optional
        x parameter for the ellipsoildal Leaf Angle Distribution function of 
        Campbell 1988 [default=1, spherical LIDF].
        叶倾角分布椭球函数的x参数
    
    Returns
    -------
    Rn_sw_veg : float
        Canopy net shortwave radiation (W m-2).
        冠层净短波辐射
    Rn_sw_soil : float
        Soil net shortwave radiation (W m-2).
        土壤净短波辐射

    References
    ----------    
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    .. [Kustas1999] Kustas and Norman (1999) Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29, http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''
    
    from math import radians, cos, sin, tan, log, sqrt, exp
    #calculate aborprtivity     计算吸收率（？）
    ameanv = 1.0-rho_leaf_vis-tau_leaf_vis
    ameann = 1.0-rho_leaf_nir-tau_leaf_nir
    # Calculate canopy beam extinction coefficient      计算冠层辐射消光系数
    #Modification to include other LADs
    akb=sqrt(x_LAD**2+tan(radians(sza))**2)/(x_LAD+1.774*(x_LAD+1.182)**-0.733) # Eq. 15.4
    # Integrate to get the diffuse transmitance     积分以获得漫射透射率

    taud=0
    for angle in range(0,90,5):
        akd=CalcKbe_Campbell(angle,x_LAD) # Eq. 15.4
        taub=exp(-akd*LAI)
        taud = taud+taub*cos(radians(angle))*sin(radians(angle))*radians(5)
    taud=2.0*taud
    #D I F F U S E   C O M P O N E N T S
    #Diffuse light canopy reflection coefficients  for a deep canopy	深冠层的漫射光冠层反射系数
    #akd=-0.0683*log(LAI)+0.804                  # Fit to Fig 15.4 for x=1
    akd=-log(taud)/LAI
    rcpyn=(1.0-sqrt(ameann))/(1.0+sqrt(ameann)) # Eq 15.7   
    rcpyv=(1.0-sqrt(ameanv))/(1.0+sqrt(ameanv))
    rdcpyn=2.0*akd*rcpyn/(akd+1.0)              #Eq 15.8      
    rdcpyv=2.0*akd*rcpyv/(akd+1.0) 
    #Diffuse canopy transmission coeff (visible)    漫射冠层透射系数(可见光)
    expfac = sqrt(ameanv)*akd*LAI
    xnum = (rdcpyv*rdcpyv-1.0)*exp(-expfac)
    xden = (rdcpyv*rsoilv-1.0)+rdcpyv*(rdcpyv-rsoilv)*exp(-2.0*expfac)
    taudv = xnum/xden                           #Eq 15.11
    #Diffuse canopy transmission coeff (NIR) 	漫射冠层透射系数(近红外)
			
    expfac = sqrt(ameann)*akd*LAI;
    xnum = (rdcpyn*rdcpyn-1.0)*exp(-expfac);
    xden = (rdcpyn*rsoiln-1.0)+rdcpyn*(rdcpyn-rsoiln)*exp(-2.0*expfac)
    taudn = xnum/xden                           #Eq 15.11
    #Diffuse radiation surface albedo for a generic canopy一般冠层的漫射面反照率

    fact=((rdcpyn-rsoiln)/(rdcpyn*rsoiln-1.0))*exp(-2.0*sqrt(ameann)*akd*LAI)   #Eq 15.9
    albdn=(rdcpyn+fact)/(1.0+rdcpyn*fact)
    fact=((rdcpyv-rsoilv)/(rdcpyv*rsoilv-1.0))*exp(-2.0*sqrt(ameanv)*akd*LAI)   #Eq 15.9
    albdv=(rdcpyv+fact)/(1.0+rdcpyv*fact)
    #B E A M   C O M P O N E N T S
    #Direct beam extinction coeff (spher. LAD)  直接光束消光系数
#==============================================================================
#     akb=0.5/coszen
#==============================================================================
    #Direct beam canopy reflection coefficients for a deep canopy   深冠层的直束冠层反射系数
    rcpyn=(1.0-sqrt(ameann))/(1.0+sqrt(ameann))                                 #Eq 15.7   
    rcpyv=(1.0-sqrt(ameanv))/(1.0+sqrt(ameanv))
    rbcpyn=2.0*akb*rcpyn/(akb+1.0)                                              #Eq 15.8      
    rbcpyv=2.0*akb*rcpyv/(akb+1.0); 
    fact=((rbcpyn-rsoiln)/(rbcpyn*rsoiln-1.0))*exp(-2.0*sqrt(ameann)*akb*LAI)   #Eq 15.9
    albbn=(rbcpyn+fact)/(1.0+rbcpyn*fact)
    fact=((rbcpyv-rsoilv)/(rbcpyv*rsoilv-1.0))*exp(-2.0*sqrt(ameanv)*akb*LAI)   #Eq 15.9
    albbv=(rbcpyv+fact)/(1.0+rbcpyv*fact)
    #Weighted average albedo 加权平均反照率
    albedo_dir=fvis*albbv+fnir*albbn
    albedo_dif=fvis*albdv+fnir*albdn
    #Direct beam+scattered canopy transmission coeff (visible)  直接光束+散射冠层传输系数(可见光)
    expfac = sqrt(ameanv)*akb*LAI
    xnum = (rbcpyv*rbcpyv-1.0)*exp(-expfac)
    xden = (rbcpyv*rsoilv-1.0)+rbcpyv*(rbcpyv-rsoilv)*exp(-2.0*expfac)
    taubtv = xnum/xden                                                          #Eq 15.11
    #Direct beam+scattered canopy transmission coeff (NIR)  直接光束+散射冠层传输系数(近红外)

    expfac = sqrt(ameann)*akb*LAI
    xnum = (rbcpyn*rbcpyn-1.0)*exp(-expfac)
    xden = (rbcpyn*rsoiln-1.0)+rbcpyn*(rbcpyn-rsoiln)*exp(-2.0*expfac)
    taubtn = xnum/xden                                                          #Eq 15.11
    tau_dir=fvis*taubtv+fnir*taubtn
    tau_dif=fvis*taudv+fnir*taudn
    albedosoil=fvis*rsoilv+fnir*rsoiln
    Rn_sw_veg=(1.0-tau_dir)*(1.0-albedo_dir)*Sdn_dir+(1.0-tau_dif)*(1.0-albedo_dif)*Sdn_dif
    Rn_sw_soil=tau_dir*(1.0-albedosoil)*Sdn_dir+tau_dif*(1.0-albedosoil)*Sdn_dif
    return Rn_sw_veg, Rn_sw_soil
