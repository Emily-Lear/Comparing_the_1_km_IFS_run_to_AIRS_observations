%% find ERA5 data closest in time to each AIRS granules in region (3 hourly) 
% and remove data outside of min/max lon/lat ranges for each AIRS granule

% Functions written by others included in this function: 
%%% getnet: Neil Hindley. University of Bath, UK. (nh351@bath.ac.uk)


% IFS hours 
ifs_h=0:3:21;
ifs_hours=permute(ifs_h,[2,1]);
ifs_hours=num2str(ifs_hours,'%.2d');

% directory for ERA5
d_era5='D:\ERA5_2018';


find_waves=load('find_waves_box2.mat');

% array of hours (IFS times) for 1st 14 days of November 2018
% (replace 13 with 9 for times in 1st 10 days of November)
hours=[3:3:21 repmat(0:3:21,1,13)]; %9
no_hours=[7 repmat(8,1,13)]; %9

Days=[];

% make day array
for Day=1:14

    hour_num=no_hours(Day);

    Days=[Days repmat(Day,1,hour_num)];

end

d_num=[];
% find MATLAB times for files 
for i=1:length(hours)

    D=Days(i);
    H=hours(i);
    
    d_num=[d_num datenum(2018,11,D,H,0,0)];

end

ifs_time=[];
era5_day=[];
% loop through AIRS granules and find closest 3 hourly ERA5 time and day for 
% ERA5
for i=1:length(find_waves.GN)

    time_m=find_waves.time_m(i);

    % find index of closest time in IFS
    [~,time_i]=min(abs(d_num-time_m));
    ifs_time=[ifs_time d_num(time_i)];

    % convert to days and hours - IFS time 
    T_i=datestr(d_num(time_i),'dd-mm HH:MM:SS');
    era5_day=[era5_day str2double(T_i(1:2))];

end


for Day=1:14 %10 (repalce 14 with 10 for 1st 10 days)

    %  find day of year 
    date=datetime(2018, 11, Day);
    d=day(date,'dayofyear');


    % ERA5 file for day
    filename=strcat('era5_2018d',num2str(d),'.nc');
    era5=fullfile(d_era5, filename);

    NC = getnet(era5);
    e_time=NC.Data.time;
    
    % Convert time to serial date number 
    era5_d=double(e_time)/24 + datenum('1900-01-01 00:00:00');

    
    % convert to days and hours - ERA5 time 
    T_era5=datestr(era5_d,'dd-mm HH:MM:SS');
    h_era5=T_era5(:,7:8);

    % find ERA5 times for IFS hours
    [~,hours_i]=ismember(str2num(ifs_hours),str2num(h_era5));
    % 3 hourly ERA5 times 
    era5_3h=era5_d(hours_i);

    % remove 1st IFS time
    if Day==1
        era5_3h(1)=[];
    end

    % find airs granule indicies closest to times in ERA5 day 
    gn_index=find(era5_day==Day);

    % loop through granules for day 
    for j=1:length(gn_index)

        % granule index
        gn_i=gn_index(j);

        GN=find_waves.GN(gn_i);
        time_m=find_waves.time_m(gn_i);
        Day_airs=find_waves.Day(gn_i);

        if Day_airs <= 14 %10 (replace with 10 for 1st 10 days)

            minLon=find_waves.minLon(gn_i);
            maxLon=find_waves.maxLon(gn_i);
            minLat=find_waves.minLat(gn_i);
            maxLat=find_waves.maxLat(gn_i);
    
            % find 3 hourly time in ERA5
            [~,t_i]=min(abs(era5_3h-ifs_time(gn_i)));
            t_3h=era5_3h(t_i);
            % find index in ERA5 times 
            time_index=find(era5_d==t_3h);
        
        
            % get ERA5 temperatures for time 
            t=squeeze(NC.Data.t(:,:,:,time_index));
        
            % ERA5 lon/lat values 
            lon=NC.Data.longitude;
            lat=NC.Data.latitude;
        
            % find lon/lat indices 
            loni=find(lon>=(minLon-2) & lon<=(maxLon+2));
            lati=find(lat>=(minLat-2) & lat<=(maxLat+2));
        
            % find area of ERA5 data 
            t=t(loni, lati,:);
        
            % lon and lat
            lon2=lon(loni);
            lat2=lat(lati);
        
            % convert to days and hours - AIRS granule time 
            T_a=datestr(time_m,'dd-mm HH:MM:SS');
            d_a=T_a(1:2);
            h_a=T_a(7:8);
            m_a=T_a(10:11);
            S_a=T_a(13:14);
        
            % convert to days and hours - ERA5 time 
            T_e=datestr(era5_d(time_index),'dd-mm HH:MM:SS');
            h_e=T_e(7:8);
            m_e=T_e(10:11);
        
            
            % save in structure 
            era5_airs.t=t;
            era5_airs.lon=lon2;
            era5_airs.lat=lat2;
            era5_airs.airs_time=time_m;
            era5_airs.level=NC.Data.level;
            era5_airs.era5_day=Day;
            era5_airs.airs_day=Day_airs;
            era5_airs.GN=GN;
            era5_airs.era5_time=era5_d(time_index);
        
            % save structure
            save(strcat(['C:\Users\ejl45\OneDrive - University of Bath\ERA5\ERA5_AIRS_3h\era5_',num2str(Day_airs),'_',h_a,'h',m_a,'m_',num2str(GN,'%.2d'),'.mat']),'-struct', 'era5_airs')
            
        end
    end

end



