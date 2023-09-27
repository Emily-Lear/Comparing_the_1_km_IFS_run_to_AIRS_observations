%% Find IFS data for AIRS granule areas and closest to AIRS times 

% Functions written by others included in this function: 
%%% getnet: Neil Hindley. University of Bath, UK. (nh351@bath.ac.uk)


ts=180:180:21600;

% find index 
i=7+9*8;
ts(i)

% file suffixes for 1 - 10th Nov 2018
times=180:180:14220;

% array of hours for 1st 14 days of November 
hours=[3:3:21 repmat(0:3:21,1,13)]; %9
no_hours=[7 repmat(8,1,13)]; %9

Days=[];

% make day array
for Day=1:14 %10 (replace with 10 for 1st 10 days of November)

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

%%

% IFS directory 
Directory='/data3/emily/1km_data';
DirectoryInfo=dir(Directory);

% filenames 
filenames=char(DirectoryInfo.name); 
filenames(1:2,:)=[];

% AIRS
find_waves=load('find_waves_box2.mat');

length(find_waves.GN);

% loop through granules 
for j=1:length(find_waves.GN)

    disp(j)
    GN=find_waves.GN(j)
    Day=find_waves.Day(j)
    time_m=find_waves.time_m(j);

    % 1st 14 days
    if Day <= 14

        % find closest time 
        % find index of closest time in IFS
        [~,time_i]=min(abs(d_num-time_m));

        % find IFS file 
        file=filenames(time_i,:);
        f=fullfile(Directory,file);

        % get IFS temperatures 
        ifs=getnet(f);

        % find region of AIRS granule 
        t=ifs.Data.t;

        % ERA5 lon/lat values 
        lon=ifs.Data.longitude;
        lat=ifs.Data.latitude;

        minLon=find_waves.minLon(j);
        maxLon=find_waves.maxLon(j);
        minLat=find_waves.minLat(j);
        maxLat=find_waves.maxLat(j);

        % find lon/lat indices 
        loni=find(lon>=(minLon-2) & lon<=(maxLon+2));
        lati=find(lat>=(minLat-2) & lat<=(maxLat+2));

        % find area of ERA5 data 
        t=t(loni, lati,:);

        % lon and lat
        lon2=lon(loni);
        lat2=lat(lati);

        % convert to days and hours - AIRS granule time 
        T_a=datestr(time_m,'dd-mm HH:MM:SS')
        d_a=T_a(1:2);
        h_a=T_a(7:8);
        m_a=T_a(10:11);
        S_a=T_a(13:14);

        % convert to days and hours - IFS time 
        T_i=datestr(d_num(time_i),'dd-mm HH:MM:SS');
        d_i=T_i(1:2);
        h_i=T_i(7:8);
        m_i=T_i(10:11);

        % save in structure 
        ifs_airs.t=t;
        ifs_airs.lon=lon2;
        ifs_airs.lat=lat2;
        ifs_airs.airs_time=time_m;
        ifs_airs.level=ifs.Data.level;
        ifs_airs.Day=Day;
        ifs_airs.ifs_Day=str2num(d_i); 
        ifs_airs.GN=GN;
        ifs_airs.ifs_time=d_num(time_i);
        

        % save structure 
        save(strcat(['/data3/emily/IFS_AIRS/ifs_',num2str(d_a),'_',num2str(h_a),'h',num2str(m_a),'m_',num2str(GN),'.mat']),'-struct', 'ifs_airs')
    
    end

end


