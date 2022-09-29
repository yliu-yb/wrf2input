from netCDF4 import Dataset
from wrf import getvar, ll_to_xy, disable_xarray,  ALL_TIMES, interplevel
import datetime
import numpy as np
from common import temperature_from_potential_temperature, calc_density

class wrftoforce():
    def __init__(self, file, zlevels, outfile):
        self.file = file
        self.zlevels = zlevels
        self.outfile = outfile
        self.DENSITY = []
        self.dt = []
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

        times = ncfile.variables["Times"]
        dt = []
        for time in times:
            strDatatime = ''
            for b in time:
                strDatatime += bytes.decode(b)
            dt.append(datetime.datetime.strptime(strDatatime, '%Y-%m-%d_%H:%M:%S'))
        self.dt = dt
        self.lon = lons[0, y, :]
        self.lat = lats[0, :, x]

        self.U = np.array(interplevel(U[:, :, :, 1:], z[:, :, :, :], self.zlevels))
        self.V = np.array(interplevel(V[:, :, 1:, :], z[:, :, :, :], self.zlevels))
        self.W = np.array(interplevel(W[:, 1:, :, :], z[:, :, :, :], self.zlevels))
        self.TH = np.array(interplevel(T[:, :, :, :], z[:, :, :, :], self.zlevels))
        self.Q = np.array(interplevel(Q[:, :, :, :], z[:, :, :, :], self.zlevels))
        self.P = np.array(interplevel(p[:, :, :, :], z[:, :, :, :], self.zlevels))

        self.DENSITY = calc_density(temperature_from_potential_temperature(self.TH, self.P), self.P * 100)

        ncfile.close()
    def write2nc(self):
        ncfile = Dataset(self.outfile, mode='w', format='NETCDF4_CLASSIC')

        ncfile.title = 'force data for Small-scale Numerical Forcast Model'
        ncfile.subtitle = self.outfile

        time_dim = ncfile.createDimension('time', None)
        lev_dim = ncfile.createDimension('height', len(self.zlevels))
        lat_dim = ncfile.createDimension('lat', len(self.lat))
        lon_dim = ncfile.createDimension('lon', len(self.lon))

        # Define two variables with the same names as dimensions,
        # a conventional way to define "coordinate variables".

        time = ncfile.createVariable('time', np.float64, ('time',))
        time.units = 's'
        time.long_name = 'seconds since 1970-01-01 00:00:00'

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
        DENSITY_bottom = ncfile.createVariable('density_bottom', np.float64,
                                          ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        DENSITY_bottom.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_bottom.standard_name = 'bottom boundary density'  # this is a CF standard name

        DENSITY_top = ncfile.createVariable('density_top', np.float64,
                                       ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        DENSITY_top.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_top.standard_name = 'top boundary density'  # this is a CF standard name

        DENSITY_west = ncfile.createVariable('density_west', np.float64,
                                        ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        DENSITY_west.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_west.standard_name = 'west boundary density'  # this is a CF standard name

        DENSITY_east = ncfile.createVariable('density_east', np.float64,
                                        ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        DENSITY_east.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_east.standard_name = 'east boundary density'  # this is a CF standard name

        DENSITY_south = ncfile.createVariable('density_south', np.float64,
                                         ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        DENSITY_south.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_south.standard_name = 'south boundary density'  # this is a CF standard name

        DENSITY_north = ncfile.createVariable('density_north', np.float64,
                                         ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        DENSITY_north.units = 'kgm-3'  # degrees kgm-3elvin
        DENSITY_north.standard_name = 'north boundary density'  # this is a CF standard name
        # ------------

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
        W_bottom = ncfile.createVariable('W_bottom', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        W_bottom.units = 'ms-1'
        W_bottom.standard_name = 'bottom boundary vertical velocity'  # Wis is a CF standard name

        W_top = ncfile.createVariable('W_top', np.float64, ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        W_top.units = 'ms-1'  
        W_top.standard_name = 'top boundary vertical velocity'  # Wis is a CF standard name

        W_west = ncfile.createVariable('W_west', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        W_west.units = 'ms-1'
        W_west.standard_name = 'west boundary vertical velocity'  # Wis is a CF standard name

        W_east = ncfile.createVariable('W_east', np.float64, ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        W_east.units = 'ms-1'
        W_east.standard_name = 'east boundary vertical velocity'  # Wis is a CF standard name

        W_south = ncfile.createVariable('W_south', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        W_south.units = 'ms-1'
        W_south.standard_name = 'south boundary vertical velocity'  # Wis is a CF standard name

        W_north = ncfile.createVariable('W_north', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        W_north.units = 'ms-1'
        W_north.standard_name = 'north boundary vertical velocity'  # Wis is a CF standard name

        # ---------------
        P_bottom = ncfile.createVariable('P_bottom', np.float64,
                                         ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        P_bottom.units = 'hPa'
        P_bottom.standard_name = 'bottom boundary pressure'  # Pis is a CF standard name

        P_top = ncfile.createVariable('P_top', np.float64,
                                      ('time', 'lat', 'lon'))  # note: unlimited dimension is leftmost
        P_top.units = 'hPa'
        P_top.standard_name = 'top boundary pressure'  # Pis is a CF standard name

        P_west = ncfile.createVariable('P_west', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        P_west.units = 'hPa'
        P_west.standard_name = 'west boundary pressure'  # Pis is a CF standard name

        P_east = ncfile.createVariable('P_east', np.float64,
                                       ('time', 'height', 'lat'))  # note: unlimited dimension is leftmost
        P_east.units = 'hPa'
        P_east.standard_name = 'east boundary pressure'  # Pis is a CF standard name

        P_south = ncfile.createVariable('P_south', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        P_south.units = 'hPa'
        P_south.standard_name = 'south boundary pressure'  # Pis is a CF standard name

        P_north = ncfile.createVariable('P_north', np.float64,
                                        ('time', 'height', 'lon'))  # note: unlimited dimension is leftmost
        P_north.units = 'hPa'
        P_north.standard_name = 'north boundary pressure'  # Pis is a CF standard name
        # -------------

        # Write latitudes, longitudes, height, time
        lat = self.lat
        lon = self.lon
        height = self.zlevels
        time[:] = [(i - datetime.datetime(1970,1,1)).total_seconds() for i in self.dt]

        # write 3d data
        DENSITY_bottom = self.DENSITY[:, 0, :, :]
        DENSITY_top = self.DENSITY[:, len(self.zlevels) - 1, :, :]
        DENSITY_west = self.DENSITY[:, :, :, 0]
        DENSITY_east = self.DENSITY[:, :, :, len(self.lat) - 1]
        DENSITY_south = self.DENSITY[:, :, 0, :]
        DENSITY_north = self.DENSITY[:, :, len(self.lon) - 1, :]

        TH_bottom = self.TH[:, 0, :, :]
        TH_top = self.TH[:, len(self.zlevels) - 1, :, :]
        TH_west = self.TH[:, :, :, 0]
        TH_east = self.TH[:, :, :, len(self.lat) - 1]
        TH_south = self.TH[:, :, 0, :]
        TH_north = self.TH[:, :, len(self.lon) - 1, :]

        U_bottom = self.U[:, 0, :, :]
        U_top = self.U[:, len(self.zlevels) - 1, :, :]
        U_west = self.U[:, :, :, 0]
        U_east = self.U[:, :, :, len(self.lat) - 1]
        U_south = self.U[:, :, 0, :]
        U_north = self.U[:, :, len(self.lon) - 1, :]
        
        V_bottom = self.V[:, 0, :, :]
        V_top = self.V[:, len(self.zlevels) - 1, :, :]
        V_west = self.V[:, :, :, 0]
        V_east = self.V[:, :, :, len(self.lat) - 1]
        V_south = self.V[:, :, 0, :]
        V_north = self.V[:, :, len(self.lon) - 1, :]
        
        W_bottom = self.W[:, 0, :, :]
        W_top = self.W[:, len(self.zlevels) - 1, :, :]
        W_west = self.W[:, :, :, 0]
        W_east = self.W[:, :, :, len(self.lat) - 1]
        W_south = self.W[:, :, 0, :]
        W_north = self.W[:, :, len(self.lon) - 1, :]

        P_bottom = self.P[:, 0, :, :]
        P_top = self.P[:, len(self.zlevels) - 1, :, :]
        P_west = self.P[:, :, :, 0]
        P_east = self.P[:, :, :, len(self.lat) - 1]
        P_south = self.P[:, :, 0, :]
        P_north = self.P[:, :, len(self.lon) - 1, :]

        # read data back from variable (by slicing it), print min and max
        print("temp.shape:", TH_bottom.shape)
        print("-- Min/Max values:", TH_bottom[:, :, :].min(), TH_bottom[:, :, :].max())
        print("-- Min/Max values:", TH_top[:, :, :].min(), TH_top[:, :, :].max())

        ncfile.close()
        pass

if __name__ == '__main__':
    zlevels = []
    zlevels.append(250)
    for i in range(0, 22):
        zlevels.append(300 + i * 100)
    zlevels.append(2450)
    wrf2force = wrftoforce('/home/yl/wrf/wrfv4.4/WRF/test/data_save/model_drive_data/wrfout_d03_2016-06-01_05:00:00', zlevels, 'force_2016-06-01_05:00:00.nc')
