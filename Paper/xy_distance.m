function [x,y] = xy_distance(Lon,Lat,Loc1)

%  Find x and y coordinates from a point using lon/lat values 

% Inputs 

% Lon: Longitude array 
% Lat: Latitude array 
% Loc1: lon/lat coordinates of point x and y distance will be found from
% e.g. [lat1 lon1]

% outputs:

% x: distance in x direction (km)
% y: distance in y direction (km) 

% Function included in this function: 
%%% nph_haversine, Neil Hindley. University of Bath, UK. 
% (nh351@bath.ac.uk)


% calculate distance and bearing 
    
% distance
lat1=Loc1(1);
lon1=Loc1(2);

a=Lat(:);
b=Lon(:);
Loc2 = double([a b]);
d=nph_haversine(Loc1,Loc2);
size(d);
d=reshape(d,size(Lon));
size(d);
dlon= Lon-lon1;

% Bearing
X =cosd(Lat).*sind(dlon);
Y= cosd(lat1).*sind(Lat)-sind(lat1).*cosd(Lat).*cosd(dlon);
       
b= atan2d(X,Y);
size(b);

% Find x and y coordinates
x=d.*sind(b); 
y=d.*cosd(b);

x=double(x);
y=double(y);

end