# sanctuary heat time-series
#
import numpy as np
import xarray as xr
from set_depth import set_depth
from shapely.geometry import Polygon, Point
import scipy.io
import seawater as sw
from matplotlib.path import Path
import pdb

url='http://oceanmodeling.pmc.ucsc.edu:8080/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd'
model=xr.open_dataset(url)

Cs_r=model['Cs_r']
Cs_r=Cs_r.values

Cs_w=model['Cs_w']
Cs_w=Cs_w.values

Vstretching=model['Vstretching']
Vtransform=model['Vtransform']

h=model['h']
hc=model['hc']

theta_b=model['theta_b']
theta_b=theta_b.values
theta_s=model['theta_s']
theta_s=theta_s.values

zeta=model['zeta']
N=42
igrid=1
report=0
# okay we want to loop through the times
[nt,nx,ny]=zeta.shape
temp=model['temp']
salt=model['salt']
# load the sanctuary data
file='c:\mbns_file\mbns_outline.mat'
mbns=scipy.io.loadmat(file)
lat=mbns['lat']
long=mbns['long']
lat=np.transpose(lat)
long=np.transpose(long)
# get the model latitude and longitude of the grid
modellat=model['lat_rho']
modellon=model['lon_rho']
#
apoly=[]
lat=lat.flatten()
long=long.flatten()
lat=np.array(lat)
long=np.array(long)
for i in range(0,6666):
    apoly.append((long[i],lat[i]))
p=Path(apoly)
flat=modellat.values.flatten()
flon=modellon.values.flatten()
points=np.vstack((flon,flat)).T
thegrid=p.contains_points(points)
sanctuarylats=flat[thegrid]
sanctuarylons=flon[thegrid]
#
nx=np.arange(0,181)
ny=np.arange(0,186)
[mx,my]=np.meshgrid(ny,nx)

mx=mx.flatten()
my=my.flatten()

xindex=mx[thegrid]
yindex=my[thegrid]
zeta=model['zeta']
#pdb.set_trace()

for k in np.arange(0,nt):
    print(k)
    azeta=np.squeeze(zeta[k,:,:])
    z=set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,igrid,h,azeta,report)
    atemp=temp[k,:,:,:]
    asalt=salt[k,:,:,:]
    atemp=atemp.values
    asalt=asalt.values
    atemp=atemp.transpose((1,2,0))
    asalt=asalt.transpose((1,2,0))    
    ztest=z[yindex,xindex,:]
    ttest=atemp[yindex,xindex,:]
    stest=asalt[yindex,xindex,:]

    heatcapacity=sw.cp(stest,ttest,ztest) # compute heat capacity of seawater
    density=sw.dens(stest,ttest,ztest) # compute density of seawhater
    # 
    [num,nz]=ztest.shape
    for j in np.arange(0,num):
        asum=0
        for i in np.arange(0,nz-1):
        # surface area is 10*10 km or 100 km^2 (ignoring the coastline)
            if ztest[j,i] < 100:
                dz=ztest[j,i+1]-ztest[j,i]
                dT=ttest[j,i+1]-ttest[j,i]
                triang=0.5*dz*dT
                partial=heatcapacity[j,i]*density[j,i]*(ttest[j,i]*dz+triang)
                asum=asum+partial
        if j==0:
            heat=asum
        else:
            heat=np.append(heat,asum)
    #
    if k==0:
        total_heat=np.sum(heat)
    else:
        total_heat=np.append(total_heat,np.sum(heat))
    #

   
