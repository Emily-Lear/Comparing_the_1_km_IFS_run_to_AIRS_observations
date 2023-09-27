function [Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing)

% Plots topography or coastlines on a grid regular in distance

% Function included in this MATLAB script:
%%% topo_v2: Corwin Wright, University of Bath UK, (c.wright@bath.ac.uk)

% Inputs: 
   % LonBox: two-element array [min lon, max lon], in range -360 to +360
   % LatBox: two-element array [min lat, max lat], in range -90 to +90
   % xRange: two-element array [min distance, max distance] (from Loc1 in x 
   %         direction)
   % yRange: two-element array [min distance, max distance] (from Loc1 in y
   %         direction)
   % Loc1:   two-element array [lat1 lon1], the coordinates from which the
   %         distance of all the other points is calculated 
   % point_spacing: Spacing between the points of the regular distance grid
   % Topo:    Topo or Coasts 

% Outputs: 
    % Topo_km:   struct containing grid of lons, lat, distance in x and y and 
    %         elevation values for the chosen region (regular in distance) 
    % Coasts: shapefile data of coastlines in the selected region
    %         (regular in distance) 


% get topography/coastlines/Image in lon/lat range 
[Topo,~,~] = topo_v2(LonBox,LatBox);

%% convert topography to regular distance grid 


%point_spacing=[15 15 1 1];

% Make regular distance grid 
xq=xRange(1):point_spacing(1):xRange(2);
yq= yRange(1):point_spacing(2): yRange(2);
[Xq,Yq]=ndgrid(xq,yq);

% Find distance from 0 
d1=sqrt(Xq.^2+Yq.^2);
arclen=km2deg(d1);

% Find azimuth
az=atan2d(Xq,Yq);

lon1=Loc1(2);
lat1=Loc1(1);

% Find lon/lat of the grid points 
[latout, lonout]=reckon(lat1,lon1,arclen,az);


%reverse order of latitudesfor topography 
[~, I]=sort(Topo.lats(:,1));
Topo.lats(I,:);

% 3D Gridded Interpolant 

%Make 3D grid 
[Xq2,Yq2]=ndgrid(xq,yq);

gridded_Topo = nan(size(Xq2));

tslice = Topo.elev(:,:);

F=griddedInterpolant(Topo.lons',Topo.lats',tslice', 'linear','none');

gridded_Topo(:,:) = F(lonout,latout);


gridded_Topo=permute(gridded_Topo,[2,1,3]);
size(gridded_Topo);

% make structure 
Topo_km=struct('elev', gridded_Topo, 'lons', lonout, 'lats', latout, 'x', Xq2, 'y',Yq2);


 

end