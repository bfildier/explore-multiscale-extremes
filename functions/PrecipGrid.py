# from myImports import *

from settings import *

import numpy as np
import pandas as pd
import xarray as xr
import math
from datetime import datetime as dt
import pickle
import copy
import gc
import warnings

class PrecipGrid():
    
    def __init__(self, casestudy, verbose=False, verbose_steps=False, overwrite=False):
        
        # copy a few attributes from casestudy object #! BF edit
        for attr in 'name','region','rel_tab_dyam_seg','model','i_t_min','i_t_max','lat_slice','lon_slice':
            
            setattr(self,attr,getattr(casestudy,attr))

        # new grid size
        if self.region == 'tropics':
            self.n_lat = 60
            self.n_lon = 360
            
        self.df = self.rel_tab_dyam_seg #! BF edit
        self.path = DIR_DYAMOND #! BF edit
        self.path_data2d = self.path #! BF edit
        self.data_names = ['CWP', 'LHF', 'OM500', 'OM850', 'Precac', 'PSFC', 'PW', 'RH500', 'SHF', 'T2mm', 'U10m', 'V10m']
        
        # Output directory
        self.path_safeguard = os.path.join(DIR_DYAMOND_PROCESSED,self.region,'SAM','regridded') #! BF edit
        os.makedirs(self.path_safeguard,exist_ok=True) #! BF edit
        
        # Combined output file
        self.outfile = 'all.nc'
        self.outfilepath = os.path.join(self.path_safeguard,self.outfile)

        # 'prepared' means: prepare_data() has been run -- to check, test if grid file is saved
        self.overwrite = overwrite
        if self.overwrite:
            self.prepared = False
        else:
            filename = 'grid.pickle'
            filepath = os.path.join(self.path_safeguard,filename)
            self.prepared = os.path.isfile(filepath)
        
        self.verbose = verbose
        self.verbose_steps = verbose_steps
        
        # could be a function to work around the dataframe with original data
        if self.model == 'SAM':
            
            self.df = self.df.drop(casestudy.times_to_ignore).reset_index(drop=True) #! BF edit
            # self.n_first_data_incomplete = sum(1-self.df["path_dyamond"].str.split(".").apply(lambda x : int(x[1][-7:-1]) >= 30000))
            # time_stamps_str = np.array([str(list(casestudy.rel_tab_dyam_seg["path_dyamond"])[i]).split('_')[-1] for i in range(len(casestudy.rel_tab_dyam_seg)-1)])
            # time_stamps_str[time_stamps_str == 'nan'] = np.nan
            # self.time_stamps = np.array(time_stamps_str,dtype=int)
            # self.df = self.df[self.time_stamps >= 30000].reset_index(drop=True)
            
            # select only chosen times -- do NOT save as self.df: error when loading Precac at i_t-1
            df_reduced = self.df[slice(self.i_t_min,self.i_t_max)].reset_index(drop=True)
            self.index_per_days = [list(np.array(group.index)+self.i_t_min) for _, group in df_reduced.groupby(['month', 'day'])]
            self.days = [dt(year=2016,month=int(group.month.iloc[0]),day=int(group.day.iloc[0])).strftime("%m-%d") for _, group in df_reduced.groupby(['month', 'day'])]
            
            self.index_per_days.remove(self.index_per_days[-1])# last day empty
            self.n_days = len(self.index_per_days)
            
    # CREATE COPY CONSTRUCTOR TO LOAD prepared grid.pickle as self
            
    def prepare_data(self):
        
        if self.prepared:
            
            filename = 'grid.pickle'
            filepath = os.path.join(self.path_safeguard,filename)
            print(filepath)
            
            if not os.path.isfile(filepath):
                
                print("ISSUE, %s does not exist"%filepath)
                pass
            
            print("importing and copying %s"%filename)
            
            # load 
            file_handler = open(filepath,'rb')
            grid = pickle.load(file_handler)
            
            # copy
            self.__dict__ = dict(grid.__dict__).copy()
            
            del grid
        
        else:

            # open data
            for file in os.listdir(self.path_data2d):
                self.template_native_df = xr.open_dataset(os.path.join(self.path_data2d,file))
                break

            # compute
            
            if self.verbose_steps: print('-- Prepare data')

            if self.verbose_steps: print('compute coord centers')
            self.lat_centers = self.template_native_df['lat'].sel(lat=self.lat_slice).values
            self.lon_centers = self.template_native_df['lon'].sel(lon=self.lon_slice).values
            self.lat_length_on_center, self.lon_length_on_center = self.__compute_length_centers_from_coord_borders__()

            if self.verbose_steps: print('compute pixel surface')
            self.pixel_surface = self.lat_length_on_center * self.lon_length_on_center      

            self.global_area = self.pixel_surface.sum()/self.n_lat/self.n_lon #depending on the remeshed grid point surface you want computed here
            self.global_lat_area = self.global_area*self.n_lon 

            self.lat_global = [i for i in range(self.n_lat)]
            self.lon_global = [j for j in range(self.n_lon)]

            ## We start by computing the area of each latitude band
            if self.verbose_steps: print('compute lat band area')
            self.lat_area = np.sum(self.pixel_surface, axis=1)
            self.cumsum_lat_area = np.cumsum(self.lat_area)

            if self.verbose_steps: print('compute i and alpha lat')
            self.i_min, self.i_max, self.alpha_i_min, self.alpha_i_max = self.__get_i_and_alpha_lat__()

            if self.verbose_steps: print('compute area by lon')
            self.slices_i_lat = [slice(i_min+1, i_max) for i_min, i_max in zip(self.i_min[:,0], self.i_max[:,0])]        
            self.area_by_lon_and_global_lat = self._compute_area_by_lon_()
            self.cumsum_area_by_lon_and_global_lat = np.cumsum(self.area_by_lon_and_global_lat, axis = 1)

            if self.verbose_steps: print('compute j and alpha lon')
            self.j_min, self.j_max, self.alpha_j_min, self.alpha_j_max = self.__get_j_and_alpha_lon__()
            self.j_min, self.j_max =self.j_min.astype(int), self.j_max.astype(int)

            if self.verbose_steps: print('build slices j lon')
            self.slices_j_lon = self.__build_slices_j_lon__()

            if self.verbose_steps: print('compute grid surface')
            self.grid_surface = self.sum_data_from_center_to_global(self.pixel_surface)
        
            # set object to prepared
            self.prepared = True
            
            # remove template data
            self.template_native_df = None  # USE GARBAGE COLLECTOR?
        
            # save itself    
            self.removeOldFile('grid.pickle')
            self.pickleSelf('grid.pickle')
        
    def removeOldFile(self,filename):
        
        filepath = os.path.join(self.path_safeguard,filename)
        if os.path.isfile(filepath):
            os.remove(filepath)

    def pickle(self,what,filename):
        
        with open(os.path.join(self.path_safeguard,filename),'wb') as f:
            pickle.dump(what, f)
        
    def pickleSelf(self, filename):
        
        self.pickle(self,filename)
        
    def compute(self): #! BF edit

        #! MC
        # there is a weird bug making the fill all.nc to be corrupted some times...
        # working around for now consist in deleting it and then recreating with build_xarray. 
        # Although weirdly the pixel surf and global pixel surf variables do not seem to survive this process
        # Calling save_mean_precip seems to corrupt it
        if os.path.isfile(os.path.join(self.path_safeguard,'all.nc')):
            self.ds = xr.open_dataset(os.path.join(self.path_safeguard,'all.nc'))
            self.ds.close()
        else:
            self.ds = self.build_xarray()
            self.create_day_dim()
            self.ds.to_netcdf(os.path.join(self.path_safeguard,'all.nc'))
            self.ds.close()

        # print(self)
        # print(self.ds)


    #! BF new
    def saveMaxPrecip(self, day0, dayf):
        
        self.savePrecipRegridded(day0,dayf,func='max')
        
    #! BF new
    def saveMeanPrecip(self, day0, dayf):
        
        self.savePrecipRegridded(day0,dayf,func='mean')
        
    #! BF new
    def savePrecipRegridded(self,day0,dayf,func='mean'):
        """
        Save to netcdf regridded stats combined
        """
        
        key = '%s_prec'%func
        
        if key in list(self.ds.variables):
            
            print('%s already computed, skipping...'%key)
        
        else : 
            
            print('%s not computed...'%key)
            
            prec_regrid = []
            for i in range(day0,dayf+1):
                
                print('computing %s for day %d, %s...'%(key,i,self.days[i]))
                da_day_i = self.regrid_and_save_by_day(i,func=func)
                prec_regrid.append(da_day_i)

            ## concat the list of dataarrays along days dimensions
            da_prec_regrid = xr.concat(prec_regrid, dim = 'days')
            
            # print('ds =',self.ds)
            # print('xr concat =',da_prec_regrid)
            
            ## add the dataarray to the dataset
            # self.ds.assign(key=da_prec_regrid)
            self.ds[key] = da_prec_regrid    
        
        # remove file on disk
        print('remove old file on disk')
        self.removeOldFile(self.outfilepath)

        # save
        print('save...')
        self.ds.to_netcdf(self.outfilepath, format='NETCDF4', mode='w')

        #close the dataset
        self.ds.close()
        
    #! BF new
    def saveVarRegridded(self,day0,dayf,varid='Prec',func='mean'):
        """
        Save to netcdf regridded stats combined
        """
        
        key = '%s_%s'%(func,varid)
        
        if key in list(self.ds.variables):
            
            print('%s already computed, skipping...'%key)
        
        else : 
            
            print('%s not computed...'%key)
            
            var_regrid = []
            for i in range(day0,dayf+1):
                
                print('computing %s for day %d, %s...'%(key,i,self.days[i]))
                da_day_i = self.regrid_and_save_by_day(i,varid=varid,func=func)
                var_regrid.append(da_day_i)

            ## concat the list of dataarrays along days dimensions
            da_var_regrid = xr.concat(var_regrid, dim = 'days')
            
            # print('ds =',self.ds)
            # print('xr concat =',da_prec_regrid)
            
            ## add the dataarray to the dataset
            # self.ds.assign(key=da_prec_regrid)
            self.ds[key] = da_var_regrid
            
        varidoutfilepath = os.path.join(self.path_safeguard,'%s.nc'%varid)  
        
        # remove file on disk
        print('remove old file on disk')
        self.removeOldFile(varidoutfilepath)

        # save
        print('save...')
        self.ds.to_netcdf(varidoutfilepath, format='NETCDF4', mode='w')

        #close the dataset
        self.ds.close()
        
    #! BF new
    def regrid_prec_by_day(self,day,func='mean'):
        """
        Compute func(precipitation) on new grid for a given day
        Return it
        """
        
        day_prec_all = []
        idx_prev = self.index_per_days[day][0]-1
        precac_prev = self.loadPrecac(idx_prev,self.df)
        
        for idx in self.index_per_days[day]:

            # Get the data for the day
            precac_current = self.loadPrecac(idx,self.df)

            # get 30mn precipitation from difference
            prec = precac_current - precac_prev
            prec.rename('precipitation')

            # compute func(precipitation)
            prec_regrid_idx = getattr(self,'spatial_%s_data_from_center_to_global'%func)(prec)
            
            # # store(first prec computed has no value)
            # if idx !=0 : day_prec_all.append(prec_regrid_idx)
            # store(first prec computed has no value)
            day_prec_all.append(prec_regrid_idx)
            
            # update prev
            precac_prev = precac_current
            
            # free up memory
            gc.collect()
            
        del prec
        del precac_current
        del precac_prev
        gc.collect()
            
        stacked_array = np.stack(day_prec_all, axis=0)
        ## hstack the data in a dataframe of same shape than
        prec_regrid = getattr(np,func)(stacked_array, axis=0) 
        prec_regrid = np.expand_dims(prec_regrid, axis=2)
        
        return prec_regrid
        
    #! BF new
    def regrid_by_day(self,day,varid='Prec',func='mean'):
        """
        Compute func(var) on new grid for a given day
        Return it
        """
        
        def regrid_time_step(idx,varid):
            
            # Get the data for the day
            var_current = self.loadVar(idx,self.df,varid=varid)

            # compute func(var)
            if var_current is not None:
                
                var_regrid_idx = getattr(self,'spatial_%s_data_from_center_to_global'%func)(var_current)
                
                del var_current
                gc.collect()

            else:
                # create an xarray filled with nans
                var_regrid_idx = self.create_empty_array()
                print('regrid nans')
                
            return var_regrid_idx
        
        
        if varid == 'Prec': # special treatment to compute difference in Precac
            
            return self.regrid_prec_by_day(day,func=func)
        
        elif varid == 'LANDMASK': # special treatment to compute only once (fixed field)
            
            # regrid only first index
            idx = self.index_per_days[day][0]
            var_regrid = regrid_time_step(idx,varid)
            var_regrid = np.expand_dims(var_regrid, axis=2)

            return var_regrid    
            
        else:
            
            day_var_all = []
            
            for idx in self.index_per_days[day]:

                # regrid time step
                var_regrid_idx = regrid_time_step(idx,varid)
            
                # store(first prec computed has no value)
                day_var_all.append(var_regrid_idx)

            stacked_array = np.stack(day_var_all, axis=0)
            ## hstack the data in a dataframe of same shape than
            var_regrid = getattr(np,'nan%s'%func)(stacked_array, axis=0) 
            var_regrid = np.expand_dims(var_regrid, axis=2)

            return var_regrid    
    
    #! BF new
    def regrid_and_save_by_day(self,day,varid='Prec',func='mean'):
        """
        Compute func(precipitation) on new grid for a given day
        Save it as netcdf file
        """
        
        filename = 'day_'+self.days[day]+'.pkl' ## +26 meanPrecip and maxprecip computed for 40 days instead of 14
        filedir = os.path.join(self.path_safeguard,'%s_%s'%(func,varid))
        os.makedirs(filedir,exist_ok=True)
        filepath = os.path.join(filedir,filename)
        
        ## load it if it exists
        if not os.path.isfile(filepath):
            
            # compute
            var_regridded = self.regrid_by_day(day,varid=varid,func=func)
            
            # save as pickle file in the directory ${func}_Prec
            with open(filepath, 'wb') as f:
                pickle.dump(var_regridded, f)
            
        else:
            
            print('%s already exists'%filepath)
            
            # load as pickle file in the directory ${func}_Prec
            with open(filepath, 'rb') as f:
                var_regridded = pickle.load(f)
            
        da_day = xr.DataArray(var_regridded, dims=['lat_global', 'lon_global', 'days'], coords={'lat_global': self.lat_global, 'lon_global': self.lon_global, 'days': [day]})
        
        return da_day
            
    
    def __repr__(self): #! BF edit
        """Creates a printable version of the Distribution object. Only prints the 
        attribute value when its string fits is small enough."""

        out = '< CaseStudy object:\n'
        # print keys
        for k in self.__dict__.keys():
            out = out+' . %s: '%k
            if k in ['dist_chunks','chunks_to_ignore']:
                # show type
                out = out+'%s\n'%str(getattr(self,k).__class__)
            else:
                if len(str(getattr(self,k))) < 80:
                    # show value
                    out = out+'%s\n'%str(getattr(self,k))
                else:
                    # show type
                    out = out+'%s\n'%getattr(self,k).__class__
        out = out+' >'

        return out
                        
    def create_day_dim(self):
        self.n_days = len(self.index_per_days)
        if type(self.n_days) is not int:
            raise ValueError('number of days is not an integer')
        self.days_dim = self.ds.expand_dims('days', None)
        self.days_values = np.arange(self.n_days)
        self.ds['days'] = xr.DataArray(np.arange(self.n_days), dims = ['days']) 

    def build_xarray(self):
        
        # dimensions
        dims = ['lat', 'lon']
        coords = {'lat': self.lat_centers, 'lon': self.lon_centers}
        da = xr.DataArray(self.pixel_surface, dims = dims, coords = coords)
        
        self.global_pixel_surface = self.sum_data_from_center_to_global(self.pixel_surface)

        dims_global = ['lat_global', 'lon_global']
        coords_global = {'lat_global': self.lat_global, 'lon_global': self.lon_global}
        da_global = xr.DataArray(self.global_pixel_surface, dims = dims_global, coords = coords_global)
        
        ds = xr.Dataset({'pixel_surf': da, 'global_pixel_surf': da_global})
                
        return ds

    def __get_i_and_alpha_lat__(self):
        i_min, i_max = np.zeros((self.n_lat, self.n_lon)), np.zeros((self.n_lat, self.n_lon))
        alpha_i_min, alpha_i_max = np.ones((self.n_lat, self.n_lon)), np.ones((self.n_lat, self.n_lon))

        for i_lat, cum_length in enumerate(self.cumsum_lat_area):
            border_left = cum_length - self.lat_area[i_lat]
            border_right = cum_length
            
            for i in range(self.n_lat):
                cum_global_length = (i+1)*self.global_lat_area
                
                if cum_global_length > border_left and ((cum_global_length < border_right) or (math.isclose(cum_global_length, border_right))):
                    bottom_contrib = (cum_global_length - border_left)/self.lat_area[i_lat]
                    top_contrib = (border_right - cum_global_length)/self.lat_area[i_lat]           
                    # if self.verbose : print('local', i_lat, cum_length,'global',  i, cum_global_length,'borders',  border_left, border_right, 'contribs', bottom_contrib, top_contrib)
                    if i != self.n_lat-1:
                        i_min[i+1, :] = i_lat
                        alpha_i_min[i+1, :] = top_contrib if not (math.isclose(cum_global_length, border_right)) else 0
                    

                    i_max[i, :] = i_lat
                    alpha_i_max[i, :] = bottom_contrib if not (math.isclose(cum_global_length, border_right)) else 1
                    
        return i_min.astype(int), i_max.astype(int), alpha_i_min, alpha_i_max

    def __get_j_and_alpha_lon__(self):
        j_min, j_max = np.zeros((self.n_lat, self.n_lon)), np.zeros((self.n_lat, self.n_lon))
        alpha_j_min, alpha_j_max = np.ones((self.n_lat, self.n_lon)), np.ones((self.n_lat, self.n_lon))
        for i in range(self.n_lat):
            cumsum_area_by_lon = self.cumsum_area_by_lon_and_global_lat[i, :]
            for j_lon, cum_length in enumerate(cumsum_area_by_lon):
                border_left = cum_length - self.area_by_lon_and_global_lat[i, j_lon]
                border_right = cum_length
                
                for j in range(self.n_lon):
                    cum_global_length = (j+1)*self.global_area
                    
                    if cum_global_length > border_left  and ((cum_global_length) < border_right or (math.isclose(cum_global_length, border_right))):

                        left_contrib = (cum_global_length - border_left)/self.area_by_lon_and_global_lat[i, j_lon]
                        right_contrib = (border_right - cum_global_length)/self.area_by_lon_and_global_lat[i, j_lon]
                        # if self.verbose : print('local', j_lon, cum_length,'global',  j, cum_global_length,'borders',  border_left, border_right, 'contribs', left_contrib, right_contrib)
                        if j!= self.n_lon-1:
                            j_min[i, j+1] = j_lon
                            alpha_j_min[i, j+1] = right_contrib if not (math.isclose(cum_global_length, border_right)) else 0
                            
                        j_max[i, j] = j_lon
                        alpha_j_max[i, j] = left_contrib if not (math.isclose(cum_global_length, border_right)) else 1
    
        return j_min, j_max, alpha_j_min, alpha_j_max    
     
    def __build_slices_j_lon__(self):
        slices_j_lon = np.empty((self.n_lat, self.n_lon), dtype=object)
        for i in range(self.n_lat):
            for j in range(self.n_lon):
                slices_j_lon[i, j] = slice(int(self.j_min[i, j])+1, int(self.j_max[i, j])) 
        return slices_j_lon
    
    def _compute_area_by_lon_(self):
            area_by_lon = np.zeros((self.n_lat, self.lon_centers.shape[0]))
            for j_lon in range(self.lon_centers.shape[0]):
                for i, slice_i_lat in enumerate(self.slices_i_lat):
                    i_min = self.i_min[i, :]
                    i_min = self.check_all_values_same(i_min)
                    i_max = self.i_max[i, :]
                    i_max = self.check_all_values_same(i_max)
                    alpha_i_min = self.alpha_i_min[i, :]
                    alpha_i_min = self.check_all_values_same(alpha_i_min)
                    alpha_i_max = self.alpha_i_max[i, :]
                    alpha_i_max = self.check_all_values_same(alpha_i_max)
                
                    ## print i_min, i_max, alpha_i_min, alpha_i_max
                    if self.verbose : print(i, i_min, i_max, alpha_i_min, alpha_i_max)
                    bottom_sum = self.pixel_surface[i_min,j_lon]*alpha_i_min
                    if self.verbose : print(bottom_sum)
                    top_sum = self.pixel_surface[i_max,j_lon]*alpha_i_max
                    if self.verbose : print(top_sum)
                    mid_sum = np.sum(self.pixel_surface[slice_i_lat, j_lon])
                    if self.verbose : print(mid_sum)
                    area_by_lon[i, j_lon] = mid_sum+bottom_sum+top_sum
                    #print everything 
                    if False : print('i', i, 'j_lon', j_lon, 'i_min', i_min, 'i_max', i_max, 'slice_i_lat', slice_i_lat, 'alpha_i_min', alpha_i_min, 'alpha_i_max', alpha_i_max, 'bottom_sum', bottom_sum, 'top_sum', top_sum, 'mid_sum', mid_sum, 'area_by_lon', area_by_lon[i, j_lon])
            
            return area_by_lon
        
    def sum_data_from_center_to_global(self, data_on_center):
        x = data_on_center
        X = np.zeros((self.n_lat, self.n_lon))
        for i, slice_i_lat in enumerate(self.slices_i_lat):
            for j, slice_j_lon in enumerate(self.slices_j_lon[i]):
                if self.verbose : print(slice_i_lat, slice_j_lon)
                mid_sum = np.sum(x[slice_i_lat, slice_j_lon])
                bottom_sum = np.sum( x[self.i_min[i,j], slice_j_lon]*self.alpha_i_min[i,j])
                top_sum = np.sum( x[self.i_max[i,j], slice_j_lon]*self.alpha_i_max[i,j])
                left_sum = np.sum( x[slice_i_lat, self.j_min[i,j]]*self.alpha_j_min[i,j])
                right_sum = np.sum( x[slice_i_lat, self.j_max[i,j]]*self.alpha_j_max[i,j])
                bottom_left_corner = x[self.i_min[i,j], self.j_min[i,j]]*self.alpha_j_min[i,j]*self.alpha_i_min[i,j]
                bottom_right_corner = x[self.i_min[i,j], self.j_max[i,j]]*self.alpha_j_max[i,j]*self.alpha_i_min[i,j]
                top_left_corner = x[self.i_max[i,j], self.j_min[i,j]]*self.alpha_j_min[i,j]*self.alpha_i_max[i,j]
                top_right_corner = x[self.i_max[i,j], self.j_max[i,j]]*self.alpha_j_max[i,j]*self.alpha_i_max[i,j]
                X[i, j] = mid_sum+bottom_sum+top_sum+left_sum+right_sum+bottom_left_corner+bottom_right_corner+top_left_corner+top_right_corner
        return X
    
    def create_empty_array(self):
        
        return np.full(self.slices_j_lon.shape,np.nan)
    
    def spatial_mean_data_from_center_to_global(self, data_on_center):
        x = data_on_center*self.pixel_surface if type(data_on_center) == np.ndarray else data_on_center.values*self.pixel_surface
        X = np.zeros((self.n_lat, self.n_lon))
        for i, slice_i_lat in enumerate(self.slices_i_lat):
            for j, slice_j_lon in enumerate(self.slices_j_lon[i]):
                if self.verbose : print(slice_i_lat, slice_j_lon)
                mid= x[slice_i_lat, slice_j_lon].flatten()
                bottom = x[self.i_min[i,j], slice_j_lon]*self.alpha_i_min[i,j].flatten()
                top = x[self.i_max[i,j], slice_j_lon]*self.alpha_i_max[i,j].flatten()
                left = x[slice_i_lat, self.j_min[i,j]]*self.alpha_j_min[i,j].flatten()
                right = x[slice_i_lat, self.j_max[i,j]]*self.alpha_j_max[i,j].flatten()
                bottom_left_corner = x[self.i_min[i,j], self.j_min[i,j]]*self.alpha_j_min[i,j]*self.alpha_i_min[i,j]
                bottom_right_corner = x[self.i_min[i,j], self.j_max[i,j]]*self.alpha_j_max[i,j]*self.alpha_i_min[i,j]
                top_left_corner = x[self.i_max[i,j], self.j_min[i,j]]*self.alpha_j_min[i,j]*self.alpha_i_max[i,j]
                top_right_corner = x[self.i_max[i,j], self.j_max[i,j]]*self.alpha_j_max[i,j]*self.alpha_i_max[i,j]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    X[i, j] = np.nansum(np.concatenate([mid ,bottom ,top ,left ,right ,
                                            np.array([bottom_left_corner ,bottom_right_corner ,top_left_corner ,top_right_corner])]))

        return X/self.grid_surface       
    
    def spatial_max_data_from_center_to_global(self, data_on_center):
        x = data_on_center if type(data_on_center) == np.ndarray else data_on_center.values
        X = np.zeros((self.n_lat, self.n_lon))
        alpha_max = self.__build_alpha_max__()
        for i, slice_i_lat in enumerate(self.slices_i_lat):
            for j, slice_j_lon in enumerate(self.slices_j_lon[i]):
                if self.verbose : print(slice_i_lat, slice_j_lon)
                m = np.nanmax(x[slice_i_lat, slice_j_lon].flatten())
                b = np.nanmax((x[self.i_min[i,j], slice_j_lon]*alpha_max[self.i_min[i,j], slice_j_lon]).flatten())
                t = np.nanmax((x[self.i_max[i,j], slice_j_lon]*self.alpha_max[self.i_max[i,j], slice_j_lon]).flatten())
                l = np.nanmax((x[slice_i_lat, self.j_min[i,j]]*self.alpha_max[slice_i_lat, self.j_min[i,j]]).flatten())
                r = np.nanmax((x[slice_i_lat, self.j_max[i,j]]*self.alpha_max[slice_i_lat, self.j_max[i,j]]).flatten())
                blc = (x[self.i_min[i,j], self.j_min[i,j]]*self.alpha_max[self.i_min[i,j], self.j_min[i,j]]).flatten()
                btc = (x[self.i_min[i,j], self.j_max[i,j]]*self.alpha_max[self.i_min[i,j], self.j_max[i,j]]).flatten()
                tlc = (x[self.i_max[i,j], self.j_min[i,j]]*self.alpha_max[self.i_max[i,j], self.j_min[i,j]]).flatten()
                trc = (x[self.i_max[i,j], self.j_max[i,j]]*self.alpha_max[self.i_max[i,j], self.j_max[i,j]]).flatten()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    X[i, j] = np.nanmax(np.array([m, b, t, l, r, blc, btc, tlc, trc],dtype=object))
        # print(X.shape, X)
        return X     
                
    def check_all_values_same(self, arr):
        first_value = arr[0]
        for value in arr:
            if value != first_value:
                raise ValueError("Array contains different values")

        return first_value  # No error raised if all values are the same

    def loadPrec(self, i_t):
        # Load DYAMOND-SAM Precac

        precac_prev = self.loadPrecac(i_t-1,self.df)
        precac_current = self.loadPrecac(i_t,self.df)
        
        # get 30mn precipitation from difference
        prec = precac_current - precac_prev
        prec.rename('precipitation')
        
        # free up memory
        del precac_prev
        del precac_current
        gc.collect()
    
        return prec     

    def loadPrecac(self, i_t, df):
        root = df.iloc[i_t]['path_dyamond']
        if type(root) != str : print(root)
        file_precac = root+'.Precac.2D.nc'
        # load
        precac = xr.open_dataarray(os.path.join(self.path_data2d, file_precac)).load().sel(lon=self.lon_slice,lat=self.lat_slice)[0]
        return precac
    
    def loadRelTable(self, which='DYAMOND_SEG'):
        
        # relation table DYAMOND-SAM -- TOOCAN segmentation masks
        if which == 'DYAMOND_SEG':
            df = pd.read_csv('/home/mcarenso/code/stage-2023-multiscale-extremes/input/relation_2_table_UTC_dyamond_segmentation.csv')
            df.sort_values(by='UTC',ignore_index=True,inplace=True)
        return df   
    
    def loadVar(self,i_t,df,varid):
    
        # get filename
        root_DYAMOND = df.iloc[i_t]['path_dyamond']

        # Load DYAMOND data
        if varid in ['LLS','LLSU','LLSV']:
            
            file_DYAMOND = root_DYAMOND+'_%s.nc'%varid
            path = os.path.join(DIR_DYAMOND_DIAG2D,file_DYAMOND)

            if os.path.exists(path):
                var_DYAMOND = xr.open_dataarray(path).load().sel(lon=self.lon_slice,lat=self.lat_slice)[0].squeeze(dim=['z'])
            else:
                var_DYAMOND = None
            
        else:
            
            file_DYAMOND = root_DYAMOND+'.%s.2D.nc'%varid
            path = os.path.join(DIR_DYAMOND,file_DYAMOND)
        
            if os.path.exists(path):
                var_DYAMOND = xr.open_dataarray(path).load().sel(lon=self.lon_slice,lat=self.lat_slice)[0]
            else:
                var_DYAMOND = None

        return var_DYAMOND
    
    def __get_coord_border_from_centers__(self, coord_centers):
        coord_borders = list()
        coord_borders.append(np.floor(coord_centers[0]))
        for i in range(len(coord_centers)-1):
            coord_borders.append((coord_centers[i]+coord_centers[i+1])/2)
        coord_borders.append(np.ceil(coord_centers[-1]))  
        return coord_borders

    def __compute_length_centers_from_coord_borders__(self):
        lat_length = np.zeros(shape=(len(self.lat_centers), len(self.lon_centers)))
        lon_length = np.zeros(shape=(len(self.lat_centers), len(self.lon_centers)))
        
        self.lat_borders = self.__get_coord_border_from_centers__(self.lat_centers)
        self.lon_borders = self.__get_coord_border_from_centers__(self.lon_centers)
        
        for i_lat in range(len(self.lat_borders)-1):
            for j_lon in range(len(self.lon_borders)-1):
                lat1, lat2, lon1, lon2 = self.lat_borders[i_lat], self.lat_borders[i_lat+1], self.lon_borders[j_lon], self.lon_borders[j_lon+1]
                lat_length[i_lat, j_lon] = self.haversine(lat1, lon1, lat2, lon1)
                lon_length[i_lat, j_lon] = self.haversine(lat1, lon1, lat1, lon2)
        return lat_length, lon_length

    def haversine(self, lat1, lon1, lat2, lon2):
        """
        Calculate the distance between two points on the Earth (specified in decimal degrees)
        using the Haversine formula.
        """
        R = 6371  # Earth's radius in kilometers

        # Convert decimal degrees to radians
        lat1_rad = math.radians(lat1)
        lon1_rad = math.radians(lon1)
        lat2_rad = math.radians(lat2)
        lon2_rad = math.radians(lon2)

        # Haversine formula
        dlat = lat2_rad - lat1_rad
        dlon = lon2_rad - lon1_rad
        a = math.sin(dlat/2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        distance = R * c

        return distance
          
    # def imshow(self, x):
    #     plt.figure(figsize=(20,10))
    #     print(x.shape)
    #     plt.imshow(x, origin = 'lower')
    #     plt.colorbar()

    def __build_alpha_max__(self):
        self.alpha_max = np.zeros(shape=(len(self.lat_centers), len(self.lon_centers)))
        for i, slice_i_lat in enumerate(self.slices_i_lat):
            for j, slice_j_lon in enumerate(self.slices_j_lon[i]):
                self.alpha_max[slice_i_lat, slice_j_lon] = 1
                self.alpha_max[self.i_min[i,j], slice_j_lon] = 1 if self.alpha_i_min[i,j] > 0.5 else 0
                self.alpha_max[self.i_max[i,j], slice_j_lon] = 1 if self.alpha_i_max[i,j] > 0.5 else 0
                self.alpha_max[slice_i_lat, self.j_min[i,j]] = 1 if self.alpha_j_min[i,j] > 0.5 else 0
                self.alpha_max[slice_i_lat, self.j_max[i,j]] = 1 if self.alpha_j_max[i,j] > 0.5 else 0
                
                self.alpha_max[self.i_min[i,j], self.j_min[i,j]] = 1 if self.alpha_i_min[i,j] > 0.5 and self.alpha_j_min[i,j] > 0.5 else 0
                self.alpha_max[self.i_max[i,j], self.j_min[i,j]] = 1 if self.alpha_i_max[i,j] > 0.5 and self.alpha_j_min[i,j] > 0.5 else 0
                self.alpha_max[self.i_min[i,j], self.j_max[i,j]] = 1 if self.alpha_i_min[i,j] > 0.5 and self.alpha_j_max[i,j] > 0.5 else 0
                self.alpha_max[self.i_max[i,j], self.j_max[i,j]] = 1 if self.alpha_i_max[i,j] > 0.5 and self.alpha_j_max[i,j] > 0.5 else 0
                
        return self.alpha_max
    
    
#####---- OLD VERSIONS -----####
    
#     # OLD VERSION, ! BF unused
#     def saveMeanPrecip(self, day0, dayf):
#         try: 
#             self.ds['mean_prec']
#         except KeyError : 
#             print('Mean prec not saved, saving...')
#             mean_prec= []
#             for i in range(self.n_days):
#                 da_day_i = self.compute_mean_prec_by_day(i)
#                 mean_prec.append(da_day_i)
#             ## concat the list of dataarrays along days dimensions
#             da_mean_prec = xr.concat(mean_prec, dim = 'days')
#             ## add the dataarray to the dataset
#             self.ds['mean_prec']=da_mean_prec
#             if os.path.isfile(os.path.join(self.path_safeguard,'all.nc')):
#                 os.remove(os.path.join(self.path_safeguard,'all.nc'))
#             self.ds.to_netcdf(os.path.join(self.path_safeguard,'all.nc'), format='NETCDF4', mode='w')
#         else : print('Mean prec already saved, skipping...')
        
#         #close the dataset
#         self.ds.close()

#     # OLD VERSION, ! BF unused
#     def saveMaxPrecip(self, day0, dayf):
#         try: 
#             self.ds['max_prec']
#         except KeyError : print('Max prec not saved, saving...')
#         else : print('Max prec already saved, loading...')
#         max_prec= []
#         for i in range(self.n_days):
#             da_day_i = self.compute_max_prec_by_day(i)
#             max_prec.append(da_day_i)
#         ## concat the list of dataarrays along days dimensions
#         da_max_prec = xr.concat(max_prec, dim = 'days')
#         ## add the dataarray to the dataset
#         self.ds['max_prec']=da_max_prec
#         if os.path.isfile(os.path.join(self.path_safeguard,'all.nc')):
#             os.remove(os.path.join(self.path_safeguard,'all.nc'))
#         self.ds.to_netcdf(os.path.join(self.path_safeguard,'all.nc'), format='NETCDF4', mode='w')
#         self.ds.close()
     
#     # OLD VERSION, ! BF unused
#     def compute_mean_prec_by_day(self, day):
#         """
#         Compute the mean precipitation for a given day
#         """
        
#         filename = 'day_'+str(day+26)+'.pkl' ## +26 meanPrecip and maxprecip computed for 40 days instead of 14
#         filedir = os.path.join(self.path_safeguard,'mean_Prec')
#         os.makedirs(filedir,exist_ok=True)
#         filepath = os.path.join(filedir,filename)
        
#         ## load it if it exists
#         if os.path.isfile(filepath):
#             with open(filepath, 'rb') as f:
#                 mean_prec = pickle.load(f)
#         else:
#             day_mean_prec = []
#             for idx in self.index_per_days[day]:
#                 # Get the data for the day
#                 data = self.loadPrec(idx)
#                 # Compute the mean precipitation
#                 mean_prec = self.spatial_mean_data_from_center_to_global(data)
#                 ## first prec computed has no value
#                 if idx !=0 : day_mean_prec.append(mean_prec)
#             stacked_mean_array = np.stack(day_mean_prec, axis=0)
#             ## hstack the data in a dataframe of same shape than
#             mean_prec = np.mean(stacked_mean_array, axis=0) 
#             mean_prec = np.expand_dims(mean_prec, axis=2)
#             ## save as pkl file in the directory mean_Prec
#             with open(filepath, 'wb') as f:
#                 pickle.dump(mean_prec, f)
                
#         da_day = xr.DataArray(mean_prec, dims=['lat_global', 'lon_global', 'days'], coords={'lat_global': self.lat_global, 'lon_global': self.lon_global, 'days': [day]})
        
#         return da_day

#     # OLD VERSION, ! BF unused
#     def compute_max_prec_by_day(self, day):
#         """
#         Compute the mean precipitation for a given day
#         """
        
#         filename = 'max_Prec/day_'+str(day+26)+'.pkl'
#         ## load it if it exists
#         if os.path.isfile(os.path.join(self.path_safeguard,filename)):
#             max_prec = pickle.load(open(os.path.join(self.path_safeguard,filename), 'rb'))
#         else:
#             max_prec = []
#             for idx in self.index_per_days[day]:
#                 # Get the data for the day
#                 data = self.loadPrec(idx)
#                 # Compute the mean precipitation
#                 day_max_prec = self.spatial_max_data_from_center_to_global(data)
#                 ## first prec computed has no value
#                 if idx !=0 : max_prec.append(day_max_prec)
#             stacked_max_array = np.stack(max_prec, axis=0)
#             ## hstack the data in a dataframe of same shape than
#             max_prec = np.max(stacked_max_array, axis=0) 
#             print(max_prec.shape)
#             max_prec = np.expand_dims(max_prec, axis=2)
#             ## save as pkl file in the directory mean_Prec
#             with open(os.path.join(self.path_safeguard,filename), 'wb') as f:
#                 pickle.dump(max_prec, f)
                
#         da_day = xr.DataArray(max_prec, dims=['lat_global', 'lon_global', 'days'], coords={'lat_global': self.lat_global, 'lon_global': self.lon_global, 'days': [day]})
        
#         return da_day
     
#     # def __repr__(self):
#     #     print('entire map', self.global_area*self.n_lon*self.n_lat)
#     #     print('tranche de latitude : ', self.global_area*self.n_lon)
#     #     print('tranche de longitude : ', self.global_area*self.n_lat)
#     #     print('s', self.global_area)
#     #     return self.name + ' ' + self.region + ' ' + self.sim

#     # OLD VERSION, ! BF unused
#     def main(self):  
#         for i_t in range(1, len(self.df)):
#             precip = self.loadPrec(i_t, self.df)
#             da_precip_i_t = xr.DataArray(precip, dims = ['lat', 'lon'], coords = {'lat': self.template_native_df['lat'].values, 'lon': self.template_native_df['lon'].values})
