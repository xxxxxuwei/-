# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:09:53 2021

@author: Lenovo
"""

import numpy as np
from osgeo import gdal
import glob
 
 
 
# 出错：SyntaxError: (unicode error) 'unicodeescape' codec can't decode bytes in position 36-37: malformed \N character escape
# 出错原因：没有在目录位置之前加r
# glob函数参考：https://blog.csdn.net/wc781708249/article/details/78185839
list_tif = glob.glob(r'G:\Desktop\大创文件\MOD09GA\*.tif')
out_path = r'G:\Desktop\大创文件\OUT\out'

for i in range(0,len(list_tif),2):

    NIR = gdal.Open(list_tif[i+1])
    RED = gdal.Open(list_tif[i])
    # (filepath, fullname) = os.path.split(NIR)
    # (prename, suffix) = os.path.splitext(fullname)

    nir = NIR.ReadAsArray() * 0.0001
    red = RED.ReadAsArray() * 0.0001

    # https://www.cnblogs.com/beile/p/10789333.html
    try:
    #   ndvi = (nir - red) / (nir + red+0.000000000001)加一个小的数能避免被0除的问题
        ndvi = (nir - red) / (nir + red)
    # ImportError: No module named XXX，捕获该类错误。
    except ZeroDivisionError:
        print('ZeroDivisionError has be found!')
 
    # 将NAN转化为0值
    nan_index = np.isnan(ndvi)
    ndvi[nan_index] = 0
    ndvi = ndvi.astype(np.float32)
    # 将计算好的NDVI保存为GeoTiff文件
    gtiff_driver = gdal.GetDriverByName('GTiff')
    # 批量处理需要注意文件名是变量，这里截取对应原始文件的不带后缀的文件名
    out_ds = gtiff_driver.Create(out_path + 'MOD09GA_NDVI_'+ str(i) + '.tif',
                                 ndvi.shape[1], ndvi.shape[0], 1, gdal.GDT_Float32)
    #  将NDVI数据坐标投影设置为原始坐标投影
    out_ds.SetProjection(NIR.GetProjection())
    out_ds.SetGeoTransform(NIR.GetGeoTransform())
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(ndvi)
    out_band.FlushCache()
