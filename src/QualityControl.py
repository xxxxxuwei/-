# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:14:02 2021

@author: Lenovo
"""

import numpy as np
import gdal
import os
import pandas as pd

def opentif(filepath):
    """
        输入文件名，返回数组、宽度、高度、仿射矩阵信息、投影信息
    :param filepath: 文件的完整路径
    :return: im_data,im_width,im_height,im_geotrans,im_proj
    """
    dataset = gdal.Open(filepath)
    im_width = dataset.RasterXSize #栅格矩阵的列数
    im_height = dataset.RasterYSize #栅格矩阵的行数
    data = dataset.ReadAsArray(0,0,im_width,im_height)#获取数据
    # data = dataset.ReadAsArray()  # 获取数据
    im_data = np.array(data)
    print("opentif的shape")
    print(im_data.shape)

    im_geotrans = dataset.GetGeoTransform()#获取仿射矩阵信息
    im_proj = dataset.GetProjection()#获取投影信息
    return(im_data,im_width,im_height,im_geotrans,im_proj)

def savetif(dataset,path,im_width,im_height,im_geotrans,im_proj):
    """
        将数组保存为tif文件
    :param dataset: 需要保存的数组
    :param path: 需要保存出去的路径，包含文件名
    :param im_width: 数组宽度
    :param im_height: 数组宽度
    :param im_geotrans: 仿射矩阵信息
    :param im_proj: 投影信息
    :return:
    """
    print(dataset)
    driver = gdal.GetDriverByName("GTiff")
    outdataset = driver.Create(path, im_width, im_height, 1, gdal.GDT_Float32)
    print (path)
    outdataset.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
    outdataset.SetProjection(im_proj)  # 写入投影
    outdataset.GetRasterBand(1).WriteArray(dataset)
    outdataset.GetRasterBand(1).SetNoDataValue(0)
    print("yes")

if __name__ == "__main__":
    inDir = r"G:\Desktop\大创文件\MOD11A1"
    # inFile = "MOD11A1.A2018001.QC_Day.tif"
    Out_Dir = r"G:\Desktop\大创文件\OUT"
    # 获取LST质量控制文件的列表
    InList_Qc = [infile for infile in os.listdir(inDir) if infile.endswith(".QC_Day.tif")]

    for InFile in InList_Qc:
        # 获取完整路径
        in_Full_Dir = inDir + os.sep + InFile
        # 打开TIF文件，获取TIF文件的信息
        InData = opentif(in_Full_Dir)
        in_Array = InData[0]
        in_Array= np.array(in_Array,dtype = np.uint8)

        print(in_Array)
        print("   ")
        # 将十进制转回到二进制
        binary_repr_v = np.vectorize(np.binary_repr)
        new = binary_repr_v(in_Array, 8)
        print(new)
        # 6-7位，是控制LST质量的字段，‘00’代表 LST error flag  <= 1k
        Error_mask = np.char.count(new,'00',start=0,end=2) == 1
        print(Error_mask)

        # 打开LST文件，获取文件名
        # G:\2018\MODIS\MOD11A1\.GeoTif_Mosaic_10000\Tif\MOD11A1.A2018002.QC_Day.tif         质量控制文件
        # G:\2018\MODIS\MOD11A1\.GeoTif_Mosaic_10000\Tif\MOD11A1.A2018002.LST_Day_1km.tif    LST文件
        in_Full_Dir_Lst = in_Full_Dir[:-10] + "LST_Day_1km.tif"
        Lst_Array = opentif(in_Full_Dir_Lst)[0]
        # 将满足质量条件的提取出来，不满足条件的设置为0，后续设置为nodata
        Out_Lst_Array = np.where(Error_mask,Lst_Array,0)
        print(Out_Lst_Array)

        print(in_Full_Dir_Lst.split("\\")[-1])
        # 将masked后的LST保存，将 0 设置为SetNoDataValue()
        savetif(Out_Lst_Array,
                Out_Dir + os.sep + in_Full_Dir_Lst.split("\\")[-1],
                InData[1],InData[2],InData[3],InData[4])
