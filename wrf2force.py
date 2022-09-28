from netCDF4 import Dataset
from wrf import getvar, ll_to_xy, disable_xarray,  ALL_TIMES, interplevel
import datetime
import numpy as np

class wrftoforce():
    def __init__(self, file, zlevels, outfile):
        self.file = file
        self.zlevels = zlevels
        self.outfile = outfile

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

        self.lon = lons[:, y, :]
        self.lat = lats[:, :, x]

        self.U = np.array(interplevel(U[:, :, :, 1:], z[:, :, :, :], self.zlevels))
        self.V = np.array(interplevel(V[:, :, 1:, :], z[:, :, :, :], self.zlevels))
        self.W = np.array(interplevel(W[:, 1:, :, :], z[:, :, :, :], self.zlevels))
        self.TH = np.array(interplevel(T[:, :, :, :], z[:, :, :, :], self.zlevels))
        self.Q = np.array(interplevel(Q[:, :, :, :], z[:, :, :, :], self.zlevels))
        self.P = np.array(interplevel(p[:, :, :, :], z[:, :, :, :], self.zlevels))

        ncfile.close()

    def write2nc(self):
        ncfile = Dataset(self.outfile, mode='w', format='NETCDF4_CLASSIC')

        ncfile.title = 'initial data for Small-scale Numerical Forcast Model'
        ncfile.subtitle = self.outfile

        time_dim = ncfile.createDimension('time', None)
        lev_dim = ncfile.createDimension('height', len(self.zlevels))
        lat_dim = ncfile.createDimension('lat', len(self.lat))
        lon_dim = ncfile.createDimension('lon', len(self.lon))

        for dim in ncfile.dimensions.items():
            print(dim)

        # Define two variables with the same names as dimensions,
        # a conventional way to define "coordinate variables".
        time = ncfile.createVariable('time', np.int32, ('time',))
        time.units = 'seconds since 1970-01-01 00:00:00'
        time.long_name = 'time'

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
        TH_bottom = ncfile.createVariable('theta_bottom', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        TH_bottom.units = 'K'  # degrees Kelvin
        TH_bottom.standard_name = 'bottom boundary potential temperature'  # this is a CF standard name

        TH_top = ncfile.createVariable('theta_top', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        TH_top.units = 'K'  # degrees Kelvin
        TH_top.standard_name = 'top boundary potential temperature'  # this is a CF standard name

        TH_west = ncfile.createVariable('theta_west', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        TH_west.units = 'K'  # degrees Kelvin
        TH_west.standard_name = 'west boundary potential temperature'  # this is a CF standard name

        TH_east = ncfile.createVariable('theta_east', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        TH_east.units = 'K'  # degrees Kelvin
        TH_east.standard_name = 'east boundary potential temperature'  # this is a CF standard name

        TH_south = ncfile.createVariable('theta_south', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        TH_south.units = 'K'  # degrees Kelvin
        TH_south.standard_name = 'south boundary potential temperature'  # this is a CF standard name

        TH_north = ncfile.createVariable('theta_north', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        TH_north.units = 'K'  # degrees Kelvin
        TH_north.standard_name = 'north boundary potential temperature'  # this is a CF standard name

        #         -------------------
        U_bottom = ncfile.createVariable('U_bottom', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        U_bottom.units = 'ms-1'  
        U_bottom.standard_name = 'bottom boundary x-component velocity'  # Uis is a CF standard name

        U_top = ncfile.createVariable('U_top', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        U_top.units = 'ms-1'  
        U_top.standard_name = 'top boundary x-component velocity'  # Uis is a CF standard name

        U_west = ncfile.createVariable('U_west', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        U_west.units = 'ms-1'  
        U_west.standard_name = 'west boundary x-component velocity'  # Uis is a CF standard name

        U_east = ncfile.createVariable('U_east', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        U_east.units = 'ms-1'  
        U_east.standard_name = 'east boundary x-component velocity'  # Uis is a CF standard name

        U_south = ncfile.createVariable('U_south', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        U_south.units = 'ms-1'  
        U_south.standard_name = 'south boundary x-component velocity'  # Uis is a CF standard name

        U_north = ncfile.createVariable('U_north', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        U_north.units = 'ms-1'  
        U_north.standard_name = 'north boundary x-component velocity'  # Uis is a CF standard name
        
        # -------------------
        V_bottom = ncfile.createVariable('V_bottom', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        V_bottom.units = 'ms-1'  
        V_bottom.standard_name = 'bottom boundary y-component velocity'  # Vis is a CF standard name

        V_top = ncfile.createVariable('V_top', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        V_top.units = 'ms-1'  
        V_top.standard_name = 'top boundary y-component velocity'  # Vis is a CF standard name

        V_west = ncfile.createVariable('V_west', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        V_west.units = 'ms-1'  
        V_west.standard_name = 'west boundary y-component velocity'  # Vis is a CF standard name

        V_east = ncfile.createVariable('V_east', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        V_east.units = 'ms-1'  
        V_east.standard_name = 'east boundary y-component velocity'  # Vis is a CF standard name

        V_south = ncfile.createVariable('V_south', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        V_south.units = 'ms-1'  
        V_south.standard_name = 'south boundary y-component velocity'  # Vis is a CF standard name

        V_north = ncfile.createVariable('V_north', np.float64, ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        V_north.units = 'ms-1'  
        V_north.standard_name = 'north boundary y-component velocity'  # Vis is a CF standard name

        # ---
        W_bottom = ncfile.createWariable('W_bottom', np.float64,
                                         ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        W_bottom.units = 'ms-1'  
        W_bottom.standard_name = 'bottom boundary vertical velocity'  # Wis is a CF standard name

        W_top = ncfile.createWariable('W_top', np.float64,
                                      ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        W_top.units = 'ms-1'  
        W_top.standard_name = 'top boundary vertical velocity'  # Wis is a CF standard name

        W_west = ncfile.createWariable('W_west', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        W_west.units = 'ms-1'
        W_west.standard_name = 'west boundary vertical velocity'  # Wis is a CF standard name

        W_east = ncfile.createWariable('W_east', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        W_east.units = 'ms-1'
        W_east.standard_name = 'east boundary vertical velocity'  # Wis is a CF standard name

        W_south = ncfile.createWariable('W_south', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        W_south.units = 'ms-1'
        W_south.standard_name = 'south boundary vertical velocity'  # Wis is a CF standard name

        W_north = ncfile.createWariable('W_north', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        W_north.units = 'ms-1'
        W_north.standard_name = 'north boundary vertical velocity'  # Wis is a CF standard name

        # ---------------
        P_bottom = ncfile.createPariable('P_bottom', np.float64,
                                         ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        P_bottom.units = 'hPa'
        P_bottom.standard_name = 'bottom boundary pressure'  # Pis is a CF standard name

        P_top = ncfile.createPariable('P_top', np.float64,
                                      ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        P_top.units = 'hPa'
        P_top.standard_name = 'top boundary pressure'  # Pis is a CF standard name

        P_west = ncfile.createPariable('P_west', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        P_west.units = 'hPa'
        P_west.standard_name = 'west boundary pressure'  # Pis is a CF standard name

        P_east = ncfile.createPariable('P_east', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        P_east.units = 'hPa'
        P_east.standard_name = 'east boundary pressure'  # Pis is a CF standard name

        P_south = ncfile.createPariable('P_south', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        P_south.units = 'hPa'
        P_south.standard_name = 'south boundary pressure'  # Pis is a CF standard name

        P_north = ncfile.createPariable('P_north', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        P_north.units = 'hPa'
        P_north.standard_name = 'north boundary pressure'  # Pis is a CF standard name
        # -------------


        # Write latitudes, longitudes, height
        lat = self.lat
        lon = self.lon
        height = self.zlevels

        # write 3d data
        TH_bottom = self.TH[:, 0, :, :]
        TH_top = self.TH[:, len(self.zlevels) - 1, :, :]

        # TH = self.TH
        # U = self.U
        # V = self.V
        # W = self.W
        # P = self.P

        # read data back from variable (by slicing it), print min and max
        print("-- Min/Max values:", TH_bottom[:, :, :].min(), TH_bottom[:, :, :].max())
        print("-- Min/Max values:", TH_top[:, :, :].min(), TH_top[:, :, :].max())

        # print("-- Min/Max values:", U[:, :, :].min(), U[:, :, :].max())
        # print("-- Min/Max values:", V[:, :, :].min(), V[:, :, :].max())
        # print("-- Min/Max values:", W[:, :, :].min(), W[:, :, :].max())
        # print("-- Min/Max values:", P[:, :, :].min(), P[:, :, :].max())

        ncfile.close()
        pass
