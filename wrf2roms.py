from netCDF4 import Dataset, date2num
from wrf import getvar, interplevel
import sys
from os.path import basename
from datetime import timedelta, date, datetime
import numpy as np
from scipy.interpolate import griddata
import math
import os
import sys
import re

def interp( srcLons,srcLats,invar2d,dstLons,dstLats ):
    py = srcLats.flatten()
    px = srcLons.flatten()
    z = np.array(invar2d).flatten()
    z[z == 1.e+37]='nan'
    X, Y = np.meshgrid(dstLons, dstLats)
    outvar2d = griddata((px,py),z, (X, Y),method='linear',fill_value=1.e+37)
    return outvar2d


def rotate( u, v, angle, missing_value):

        # For each j (eta axis)...
        for j in range(len(v)):
        
            # For each i (xi axis)...
            for i in range(len(u[j])):
            
                # Check if all velues are not NaN and not a missing value
                if (
                        Double.isNaN(u[j][i])==false &&
                                Double.isNaN(v[j][i])==false &&
                                Double.isNaN(angle[j][i]) == false &&
                                u[j][i] != missingValue &&
                                v[j][i] != missingValue &&
                                angle[j][i] != missingValue
                ):
                    rotU = (u[j][i] * Math.cos(angle[j][i]) + v[j][i] * Math.sin(angle[j][i]));
                    rotV = (v[j][i] * Math.cos(angle[j][i]) - u[j][i] * Math.sin(angle[j][i]));

                    u[j][i]=rotU;
                    v[j][i]=rotV;
                }
            }

        }

    }


if len(sys.argv)!=4:
  print("Usage: python "+str(sys.argv[0])+" source_path grid_filename destination_file")
  sys.exit(-1)

src_path=sys.argv[1]
grid_filename=sys.argv[2]
dst=sys.argv[3]

try:
  os.remove(dst)
except:
  pass



# Create an empty destination file
ncdstfile = Dataset(dst,"w", format="NETCDF4")

u10mVar = ncdstfile.createVariable("U10M", "f4",("time","latitude","longitude"),fill_value=1.e+37)
u10mVar.description="grid rel. x-wind component"
u10mVar.standard_name="u-component"
u10mVar.units="m s-1"


# Open the NetCDF domain grid file
ncgridfile = Dataset(grid_filename)

RHOlat = np.array(getvar(ncgridfile, "lat_rho",meta=False))
RHOlon = np.array(getvar(ncgridfile, "lon_rho",meta=False))

angle = np.array(getvar(ncgridfile, "angle",meta=False))

print ("Interpolating...")

sys.exit(0)

srcs=os.listdir(src_path)
t=0
for src in srcs:

  print(t,src)

  # Open the NetCDF file
  ncsrcfile = Dataset(src_path+"/"+src)

  Xlat = np.array(getvar(ncsrcfile, "XLAT",meta=False))
  Xlon = np.array(getvar(ncsrcfile, "XLONG",meta=False))

  datetimeStr = str(b"".join(ncsrcfile.variables["Times"][:][0])).split("_")
  dateStr=str("".join(datetimeStr[0].split("-"))).replace("b'","")
  timeStr=str("".join(datetimeStr[1].split(":"))).replace("'","")


  # Get the wind at 10m u and v components (meteo oriented)
  uvmet10=getvar(ncsrcfile, "uvmet10",meta=False)

  print ("uvmet10:",uvmet10.shape)
  

  u10m = interp(Xlon,Xlat,uvmet10[0],RHOlon,RHOlat)
  v10m = interp(Xlon,Xlat,uvmet10[1],RHOlon,RHOlat)

  rotate(u10m,v10m,angle,1.e+37)

  print ("u10m:",u10m.shape)

  
  


