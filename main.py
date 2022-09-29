from wrf2initial import wrftoinitial
from wrf2force import wrftoforce
import os

wrfout_path = '/home/yl/wrf/wrfv4.4/WRF/test/data_save/model_drive_data_2'
wrfout_files = []

for root, dirs, files in os.walk(wrfout_path):
    for file in files:
        wrfout_files.append(file)
wrfout_files = sorted(wrfout_files)

zlevels = []
zlevels.append(250)

for i in range(0, 22):
    zlevels.append(300 + i * 100)
zlevels.append(2450)

for file in wrfout_files:
    print(file)
    wrftoinitial(wrfout_path + '/' + file, '' , zlevels, 'init_' + file)
    wrftoforce(wrfout_path + '/' + file, zlevels, 'force_' + file)

