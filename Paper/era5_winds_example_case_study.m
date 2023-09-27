%% u and v components of winds in ERA5 for 5th and 9th November 


f_05="D:\ERA5_winds\era5_2018d309.nc";
f_09="D:\ERA5_winds\era5_2018d313.nc";
files=[f_05 f_09];
days=["05" "09"];

%NC=getnet(f);
nci=ncinfo(f_05);

% NC_05=getnet(f_05);
% NC_09=getnet(f_09);

nc_days=["NC_05" "NC_09"];
wind_uv=["u" "v"];

%%
HeightRange=[26 55];
point_spacing=[30 30 1];
xRange=[-1700 2200];
yRange=[-2200 2200];
lat1=52;
lon1=94;

% loop through files 
for day_i=1:2

    f=files(day_i);
    NC=getnet(f);
    day=days(day_i);

    lon=NC.Data.longitude;
    lat=NC.Data.latitude;
    [Lon,Lat]=ndgrid(lon,lat);
        

    % wind component
    for wind_i=1:2

        wind1=wind_uv(wind_i);
        
        % Convert to serial date number 
        date=double(NC.Data.time)/24 + datenum('1900-01-01 00:00:00');
        Time=datestr(date,'dd-mmm-yyyy HH:MM');
        
        
        NC.(wind1)=squeeze(NC.Data.(wind1)(:,:,:,3));

        
        % figure
        % t=tiledlayout(2,2);
        % t.TileSpacing='tight';
        % t.Padding='tight';
        % set(gcf,'color','w');
        % 
        % ax(1)=nexttile(1);
        % pcolor(Lon_05, Lat_05, NC_05.u); shading flat
        % 
        % ax(2)=nexttile(2);
        % pcolor(Lon_05, Lat_05, NC_05.v); shading flat
        % 
        % ax(3)=nexttile(3);
        % pcolor(Lon_09, Lat_09, NC_09.u); shading flat
        % 
        % ax(4)=nexttile(4);
        % pcolor(Lon_09, Lat_09, NC_09.v); shading flat
        
        
        %% interpolate data to regular distance grid 
        

        
        % Altitude of levels 
        [~,Altitude] = ecmwf_prs_v3(137);
        z=Altitude;
        
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
        w=squeeze(NC.(wind1)(min_loni:max_loni,min_lati:max_lati,:));
        
        clear min_loni max_loni min_lati max_lati
        
        w1=w(:,I,:,:);
        
        clear w
        
        %% 3D Gridded Interpolant 
         
        sz=size(w1);
        
        % Make 3D grid 
        [Xq2,~,~]=ndgrid(xq,yq,z);
        
        griddedt = nan(size(Xq2));
        
        clear Xq2
        
        for z1 = 1:sz(3)
        
            wslice = w1(:,:,z1);
            
            F=griddedInterpolant(lon_grid,lat_grid,wslice, 'linear','none');
        
            griddedw(:,:,z1) = F(lonout,latout);
        
        end
        
        clear w1 lon_grid lat_grid
        
        %griddedt=permute(griddedt,[2,1,3]);
        size(griddedw);
        
        %% Interpolate Altitude 
        
        % swap dimensions
        wind=permute(griddedw,[3,1,2]);
        size(wind);
        
        clear griddedw
        
        zq2=HeightRange(1):point_spacing(3):HeightRange(2);
        %z=Altitude(7:50);
        
        % interpolate 
        grid_w=interp1(z,wind,zq2,'linear');
        size(grid_w);
        clear wind
        
        griddedw1.(wind1)=permute(grid_w,[2,3,1]);
        size(griddedw1);
        clear grid_w

    end

    griddedw1.x=xq; 
    griddedw1.y=yq;
    griddedw1.z=zq2;
    
    % save structure 
    save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\example case study ERA5 winds\',day,'_ERA5_winds.mat'),'-struct', 'griddedw1');
        
end
