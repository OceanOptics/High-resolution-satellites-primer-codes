from polymer.main import run_atm_corr
from polymer.level1 import Level1
from polymer.level2 import Level2
from polymer.level1_msi import Level1_MSI
from polymer.level2_nc import Level2_NETCDF
import os


if __name__ == "__main__":
    #baseLoc = '/Users/Gabe/Documents/GitHub/OSI/testfiles/'
    #fNames = ["LC08_L1TP_011029_20160823_20180130_01_T1", "LC08_L1TP_011029_20200615_20200625_01_T1"]
    #inFolders = ["1129", "1129"]
    #outFolders = ["20160823", "20200615"]
    #for (fileName, inFolder, outFolder) in zip(fNames, inFolders, outFolders):
    #    fileLoc = baseLoc + inFolder + '/' + fileName + '/'
    #   outLoc = baseLoc + inFolder + '/' + outFolder + '/polymer'
    run_atm_corr(
    Level1_MSI('/home/bjiang/Sentinel2/S2A_MSIL1C_20160803T153602_N0204_R111_T19TDJ_20160803T153929.SAFE/GRANULE/L1C_T19TDJ_A005828_20160803T153929' ,resolution='20'),
    Level2_NETCDF(outdir='/home/bjiang/Sentinel2', ext='20mpolymer.nc',
                  overwrite=True), multiprocessing = -1
    )




