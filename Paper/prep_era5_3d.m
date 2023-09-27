function [ERA5] = prep_era5_3d(f,lon1,lat1,HeightRange, xRange, yRange, time_i, point_spacing)

% Interpolate IFS data on a regular lon/lat grid to a regular distance grid
% and find temperature perturbations

% Functions written by others included in this function: 
%%% ecmwf_prs_v3, alt2pres_complex: Corwin Wright. University of Bath, UK.
% (cw785@bath.ac.uk)
%%% getnet: Neil Hindley. University of Bath, UK. (nh351@bath.ac.uk)
%%% smoothn: Anil Gannepali. 25/NOV/2001 (https://uk.mathworks.com/matlabcentral/fileexchange/725-smoothn)


% inputs

% f: file path for ERA5 data 
% lon1: Centre longitude for area 
% lat1: Centre latitude for area of data chosen
% HeightRange: Altitude range
% xRange: distance range in km
% yRange: distance range in km 
% time_i: index of time ('hours since 1900-01-01 00:00')
% point_spacing: spacing between points for all dimensions e.g. [50 50 1 1];

% Outputs

% ERA5: Structure 


ERA5=struct();

% Altitude of levels 
[~,Altitude] = ecmwf_prs_v3(137);
z=Altitude;

NC = getnet(f);

% latitudes and longitudes of points 
lon=NC.Data.longitude;
lat=NC.Data.latitude; 
t=NC.Data.t;

% if longitude range is 0 - 360
if max(lon)==360
    % change lon range to -180 - 180
    lon(lon>180)=lon(lon>180)-360;
end

clear NC

% Make regular distance grid 
xq=xRange(1):point_spacing(1):xRange(2);
yq= yRange(1):point_spacing(2): yRange(2);
[Xq,Yq]=ndgrid(xq,yq);

% Find distance from 0 
d1=sqrt(Xq.^2+Yq.^2);
arclen=km2deg(d1);

% Find azimuth
az=atan2d(Xq,Yq);

clear Xq Yq

% Find lon/lat of the grid points 
[latout, lonout]=reckon(lat1,lon1,arclen,az);

clear arclen az

% Find lon/lat range for regular distance grid 
min_latout=min(latout,[],'all');
max_latout=max(latout,[],'all');
min_lonout=min(lonout,[],'all');
max_lonout=max(lonout,[],'all');


% latitudes and longitudes of data points 
lati=find(lat>=min_latout-2 & lat<=max_latout+2);
loni=find(lon>=min_lonout-2 & lon<=max_lonout+2);
lat2=lat(lati);
lon2=lon(loni);
max_lati=max(lati);
min_lati=min(lati);
max_loni=max(loni);
min_loni=min(loni);

clear min_lonout max_lonout min_latout max_latout

% reverse order of latitudes 
[lat3, I]=sort(lat2);

clear lat2

% make a grid of lon/lat values
[lon_grid, lat_grid]=ndgrid(lon2, lat3);

clear lon2

% select temperature data 
t=squeeze(t(min_loni:max_loni,min_lati:max_lati,:,time_i));

clear min_loni max_loni min_lati max_lati

t1=t(:,I,:,:);

clear t

%% 3D Gridded Interpolant 
 
sz=size(t1);

% Make 3D grid 
[Xq2,~,~]=ndgrid(xq,yq,z);

griddedt = nan(size(Xq2));

clear Xq2

for z1 = 1:sz(3)

    tslice = t1(:,:,z1);
    
    F=griddedInterpolant(lon_grid,lat_grid,tslice, 'linear','none');

    griddedt(:,:,z1) = F(lonout,latout);

end

clear t1 lon_grid lat_grid

size(griddedt);

%% Interpolate Altitude 

% swap dimensions
temp=permute(griddedt,[3,1,2]);
size(temp);

clear griddedt

zq2=HeightRange(1):point_spacing(3):HeightRange(2);

% interpolate 
grid_t=interp1(z,temp,zq2,'linear');
size(grid_t);
clear temp

griddedt1=permute(grid_t,[2,3,1]);
size(griddedt1);
clear grid_t

%% Find Background

sz5=[11 11 1];
length(sz5)

bg = smoothn(griddedt1,sz5,'gaussian');

%% Find temperature perturbation 
tp=griddedt1(:,:,:)-bg;
size(tp);
sz_zq2=size(zq2);

lon3=repmat(lonout,[1 sz_zq2]);
lat3=repmat(latout,[1 sz_zq2]);
[~,~,z3]=ndgrid(xq, yq, zq2);
[x2,y2]=ndgrid(xq, yq);

sz = size(z3);
z3 = z3(:);

p = alt2pres_complex(z3);
p = reshape(p, sz);

ERA5.Tp=tp;
ERA5.lon=lon3;
ERA5.lat=lat3; 
ERA5.x= x2;
ERA5.y= y2;
ERA5.z= zq2;
ERA5.P=p;
ERA5.Tb=bg;
ERA5.T=griddedt1;


end



