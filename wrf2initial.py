from netCDF4 import Dataset
from wrf import getvar, ll_to_xy, disable_xarray,  ALL_TIMES, interplevel
import datetime
import numpy as np

class wrftoinitial():
    def __init__(self, file, dt, zlevels, outfile):
        self.file = file
        self.dt = dt
        self.outfile = outfile
        self.zlevels = zlevels
        self.U = []
        self.V = []
        self.W = []
        self.TH = []
        self.Q = []
        self.P = []
        self.lon = []
        self.lat = []
        self.getdata()
        self.write2nc()

    def getdata(self):
        ncfile = Dataset(self.file)
        disable_xarray()
        x, y = ll_to_xy(ncfile, 32.74, 117.14)
        z = getvar(ncfile, "z", timeidx=ALL_TIMES)
        p = getvar(ncfile, "pressure", timeidx=ALL_TIMES)
        lats = getvar(ncfile, "lat", timeidx=ALL_TIMES)
        lons = getvar(ncfile, "lon", timeidx=ALL_TIMES)
        U = getvar(ncfile, 'U', timeidx=ALL_TIMES)
        V = getvar(ncfile, 'V', timeidx=ALL_TIMES)
        W = getvar(ncfile, 'W', timeidx=ALL_TIMES)
        T = getvar(ncfile, 'T', timeidx=ALL_TIMES)
        Q = getvar(ncfile, 'QVAPOR', timeidx=ALL_TIMES)
        T = T + 300
        t_idx = 0
        self.lon = lons[t_idx,y,:]
        self.lat = lats[t_idx,:,x]

        self.U = np.array(interplevel(U[t_idx, :, :, 1:], z[t_idx,:,:,:], self.zlevels))
        self.V = np.array(interplevel(V[t_idx, :, 1:, :], z[t_idx,:,:,:], self.zlevels))
        self.W = np.array(interplevel(W[t_idx, 1:, :, :], z[t_idx,:,:,:], self.zlevels))
        self.TH = np.array(interplevel(T[t_idx, :, :, :], z[t_idx,:,:,:], self.zlevels))
        self.Q = np.array(interplevel(Q[t_idx, :, :, :], z[t_idx,:,:,:], self.zlevels))
        self.P = np.array(interplevel(p[t_idx, :, :, :], z[t_idx,:,:,:], self.zlevels))

        ncfile.close()
    def write2nc(self):
        ncfile = Dataset(self.outfile, mode='w', format='NETCDF4_CLASSIC')

        ncfile.title = 'initial data for Small-scale Numerical Forcast Model'
        ncfile.subtitle = self.outfile

        lev_dim = ncfile.createDimension('height', len(self.zlevels))
        lat_dim = ncfile.createDimension('lat', len(self.lat))
        lon_dim = ncfile.createDimension('lon', len(self.lon))

        # Define two variables with the same names as dimensions,
        # a conventional way to define "coordinate variables".
        lat = ncfile.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        height = ncfile.createVariable('height', np.float32, ('height',))
        height.units = 'm'
        height.long_name = 'geography height'

        # Define a 3D variable to hold the data
        TH = ncfile.createVariable('theta', np.float64, ('height', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        TH.units = 'K'  # degrees Kelvin
        TH.standard_name = 'potential temperature'  # this is a CF standard name

        U = ncfile.createVariable('U', np.float64,
                                   ('height', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        U.units = 'ms-1'  # degrees Kelvin
        U.standard_name = 'x-component wind'  # this is a CF standard name

        V = ncfile.createVariable('V', np.float64,
                                  ('height', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        V.units = 'ms-1'  # degrees Kelvin
        V.standard_name = 'y-component wind'  # this is a CF standard name

        W = ncfile.createVariable('W', np.float64,
                                  ('height', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        W.units = 'ms-1'  # degrees Kelvin
        W.standard_name = 'vertical velocity'  # this is a CF standard name

        P = ncfile.createVariable('pressure', np.float64,
                                  ('height', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        P.units = 'hPa'
        P.standard_name = 'pressure'  # this is a CF standard name

        # Write latitudes, longitudes, height
        lat = self.lat
        lon = self.lon
        height = self.zlevels

        # write 3d data
        TH = self.TH
        U = self.U
        V = self.V
        W = self.W
        P = self.P

        # read data back from variable (by slicing it), print min and max
        print("-- Min/Max values:", TH[:, :, :].min(), TH[:, :, :].max())
        print("-- Min/Max values:", U[:, :, :].min(), U[:, :, :].max())
        print("-- Min/Max values:", V[:, :, :].min(), V[:, :, :].max())
        print("-- Min/Max values:", W[:, :, :].min(), W[:, :, :].max())
        print("-- Min/Max values:", P[:, :, :].min(), P[:, :, :].max())

        ncfile.close()

if __name__ == '__main__':
    zlevels = []
    zlevels.append(250)
    for i in range(0, 22):
        zlevels.append(300 + i * 100)
    zlevels.append(2450)
    wrf2ini = wrftoinitial('/home/yl/wrf/wrfv4.4/WRF/test/data_save/model_drive_data/wrfout_d03_2016-06-01_05:00:00', '', zlevels, 'init_2016-06-01_05:00:00.nc')
    pass