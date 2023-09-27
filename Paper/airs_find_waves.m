function [waves] = airs_find_waves(Year, Month, DayRange,Directory,LatRange, LonRange, LatRange2, LonRange2)

% Function Description: 

% Finds variance of AIRS granules and mean longitude and latitude within a
% specific latitude/longitude range for AIRS and a range of
% days 
% outputs listed in order from highest to lowest variance 

% Function included in this function: 
%%% prep_airs_3d: Corwin Wright. University of Bath, UK. (cw785@bath.ac.uk)

% Inputs:

% Year: year to load AIRS granules for 
% Month: number of month to load AIRS granules for e.g.11
% DayRange: Range of days in month to load AIRS granules for
% e.g. [1:14]
% Directory: file path for AIRS granules e.g. 'file_path/AIRS/'
% LatRange: mean latitude range for granules e.g. [0 90]
% LonRange: mean longitude range for granules e.g.[0 180]
% LatRange2: range of latitudes to find max and min lon between for
% granules
% LonRange2: range of longitudes to find min and max lon in

% Outputs:

% waves: structure 

% Fields:
% Var: Variance of the AIRS granule 
% mLat: geographic mean Latitude of AIRS granule 
% mLon: geographic mean Longitude of AIRS granule 
% minLon: minumum longitude value
% maxLon: maximum lontitude value
% minLat: minumum latitude value
% maxLat: maximum " "
% GN: Granule number of Airs granule 
% Day: AIRS granule day
% time_m: Mean time 


% Make empty arrays of NaN values to fill

waves.Var=NaN(length(DayRange),240); % variance
waves.mLat=NaN(length(DayRange),240); % mean latitude 
waves.mLon=NaN(length(DayRange),240); % mean longitude 
waves.minLon=NaN(length(DayRange),240); % min longitude
waves.maxLon=NaN(length(DayRange),240); % max longitude
waves.minLat=NaN(length(DayRange),240); % min latitude
waves.maxLat=NaN(length(DayRange),240); % max latitude

waves.GN=NaN(length(DayRange),240); % Granule number 
waves.Day=NaN(length(DayRange), 240); % day of month
waves.time_m=NaN(length(DayRange), 240); % time 

% loop through days 
for i=1:length(DayRange)

    % find day of the year 
    t=datetime(Year, Month, DayRange(i));
    DayofYear=day(t,'dayofyear');
    
    % file name 
    Directory2=strcat(Directory,'/',num2str(Year),'/',num2str(DayofYear));

    % Find number of granules for each day
    d= dir(Directory2);
    %clear Directory2
    cell_d={d.name};
    clear d

    % delete 1st 2 elements 
    cell_d(1:2)=[];
    file_name=char(cell_d);

    % loop through granules 
    for j=1:length(cell_d)

        GN=str2double(file_name(j,15:17));

        % Airs file path 
        f= strcat(Directory2,'/',file_name(j,:));

        % check temperature values are not all nan 
        T=ncread(f,'ret_temp');
        if isnan(max(T(:),[],'omitnan'))

            disp([f,' All temperature values are NaN'])
        else

        % get longitudes and latitudes 
        lon=ncread(f,'l1_lon');
        lat=ncread(f,'l1_lat');

        clear f

        % Find mean geographic coordinates 
        [m_Lat, m_Lon]= meanm(lat(:), lon(:));
        

        % find if any lon/lat coordinates are within lon/lat range 
        if any(lat(:)>LatRange(1)-1 & lat(:)<LatRange(2)+1 & lon(:)>LonRange(1)-1 & lon(:)<LonRange(2)+1)

            [Airs] = prep_airs_3d(datenum(Year,Month,DayRange(i)),GN,'fulldatadir',Directory); 
        
            waves.time_m(i,j)=mean(Airs.l1_time(:),'omitnan');
  
            % find level closest to altitude in km 
            alt= 39;
            [~,level_i]=min(abs(Airs.ret_z-alt));
            Airs.Tp=Airs.Tp(:,:,level_i);
            % Find variance of AIRS granule at 39 km 
            waves.Var(i,j)=var(Airs.Tp(:));

            % find minimum and maximum latitudes and longitudes within
            % range
            lat_airs_i=find(Airs.l1_lat>=LatRange2(1) & Airs.l1_lat<=LatRange2(2));
            lon_airs_i=find(Airs.l1_lon>=LonRange2(1) & Airs.l1_lon<=LonRange2(2));
            minLat=min(Airs.l1_lat(lat_airs_i));
            maxLat=max(Airs.l1_lat(lat_airs_i));
            minLon=min(Airs.l1_lon(lon_airs_i));
            maxLon=max(Airs.l1_lon(lon_airs_i));


        
            clear Airs lon lat
        
            % Fill empty arrays 
            waves.GN(i,j)=GN;
            waves.Day(i,j)= DayRange(i);
            waves.mLon(i,j)=m_Lon;
            waves.mLat(i,j)=m_Lat;
            waves.minLat(i,j)=minLat;
            waves.maxLat(i,j)=maxLat;
            waves.minLon(i,j)=minLon;
            waves.maxLon(i,j)=maxLon;
        
            clear m_Lat m_Lon minLon maxLon minLat maxLat
        end

        end
    end
end

% Find NaN value index in flattened Variance 
nan_i=find(isnan(waves.Var(:))); 

% Flatten arrays % Remove values where vaiance is NaN
waves.Var=waves.Var(:); waves.Var(nan_i)=[]; 
waves.mLat=waves.mLat(:); waves.mLat(nan_i)=[];
waves.mLon=waves.mLon(:); waves.mLon(nan_i)=[];

waves.minLat=waves.minLat(:); waves.minLat(nan_i)=[];
waves.maxLat=waves.maxLat(:); waves.maxLat(nan_i)=[];
waves.minLon=waves.minLon(:); waves.minLon(nan_i)=[];
waves.maxLon=waves.maxLon(:); waves.maxLon(nan_i)=[];

waves.GN=waves.GN(:); waves.GN(nan_i)=[];
waves.Day=waves.Day(:); waves.Day(nan_i)=[];
waves.time_m=waves.time_m(:); waves.time_m(nan_i)=[];

clear nan_i

% Sort all variables in order of descending variance 
[waves.Var, I]=sort(waves.Var,'descend');
waves.mLat=waves.mLat(I);
waves.mLon=waves.mLon(I);

waves.minLat=waves.minLat(I);
waves.maxLat=waves.maxLat(I);
waves.minLon=waves.minLon(I);
waves.maxLon=waves.maxLon(I);

waves.GN=waves.GN(I);
waves.Day=waves.Day(I);
waves.time_m=waves.time_m(I);

clear I


end