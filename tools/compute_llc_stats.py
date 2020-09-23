import glob
import xarray as xr

from model_diagnostics import *

data_root = '/data/bia/synthetic_flow/llc4320/'
var_list = ['hs', 'dp', 'fp', 'spr', 'dir', 't0m1']

##########################
output = []
diagnostic_functions = [basic_stats]
for var in var_list:
    grid_files = glob.glob(data_root+'gridded/*%s.nc' %var)
    for f in grid_files:
        output.append(analize_llc_member(f, var, diagnostic_functions))
        print("processing %s" %os.path.basename(f))
var = 'hs'
diagnostic_functions = [hs_spectral_slope]
grid_files = glob.glob(data_root+'gridded/*%s.nc' %var)
for f in grid_files:
    output.append(analize_llc_member(f, var, diagnostic_functions))
    print("processing %s" %os.path.basename(f))
var = 'cur'
diagnostic_functions = [flow_stats]
grid_files = glob.glob(data_root+'gridded/*%s.nc' %var)
for f in grid_files:
    output.append(analize_llc_member(f, var, diagnostic_functions))
    print("processing %s" %os.path.basename(f))

ds = xr.merge(output)
df = ds.to_dataframe()
df = df.reset_index()
data = df.to_xarray()
data.to_netcdf(path='../data/model_stats/llc_gridded_stats.nc', mode='w')
