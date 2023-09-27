%% Resample IFS as AIRS (imgaussfilt3() function) -> noise
% 3x3 presmoothing

% Functions written by others included in this MATLAB script:
%%% getnet, nph_2dst_plus1 : Neil Hindley. University of Bath, UK.(nh351@bath.ac.uk)
%%% ecmwf_prs_v3, prep_airs_3d, alt2pres_complex, which_airs_retrieval, 
% cjw_airdensity: Corwin Wright. University of Bath, UK. (cw785@bath.ac.uk) 
%%% smoothn: Anil Gannepali. 25/NOV/2001 (https://uk.mathworks.com/matlabcentral/fileexchange/725-smoothn)


% Directory of IFS file for AIRS granules 
d_ifs='D:\IFS_AIRS';
DirectoryInfo=dir(d_ifs);
filenames=char(DirectoryInfo.name);
filenames(1:2,:)=[];

% AIRS files directory 
a='C:\Users\ejl45\OneDrive - University of Bath\AIRS\2018\';

% AIRS noise directory
d_airs_noise='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise new';

HeightRange=[26 55];
point_spacing=[2.7 2.7 0.1];
lat1=52;
lon1=94;

c=[0.25 0.25 0.25];
minwavelengths= [60 60 6];
maxwavelengths= [800 800 45];
nfreqs=500;

for i=1:length(filenames) 

    disp(i)

    % IFS file
    f_ifs=fullfile(d_ifs, strcat(filenames(i,:)));


    ifs=load(f_ifs);

    Day=ifs.Day;
    % if day is >=11
    if Day>=11
        disp('>=11')
    GN=ifs.GN;

    clear ifs.Day ifs.GN

    % find time in hours 
    T_ifs=datestr(ifs.ifs_time,'dd-mm HH:MM:SS');
    h_ifs=T_ifs(7:8);
    m_ifs=T_ifs(10:11);
    S_ifs=T_ifs(13:14);

    % find airs time in hours 
    T_airs=datestr(ifs.airs_time,'dd-mm HH:MM:SS');
    h_airs=T_airs(7:8);
    m_airs=T_airs(10:11);
    S_airs=T_airs(13:14);


    %  find day of year (airs) 
    date=datetime(2018, 11, Day);
    d=day(date,'dayofyear');

    % AIRS directory for day of year
    d_a=strcat(a, num2str(d));

    % add 0s in front of days and GN
    Day=num2str(Day,'%02d');
   
    GN=num2str(GN,'%03d');

    filename=strcat('airs_2018_',num2str(d),'_', GN,'.nc');
    f_a=fullfile(d_a, filename);

    airs=getnet(f_a);

    %% Put IFS onto regular distance grid 
  
    % Altitude of levels 
    [~,Altitude] = ecmwf_prs_v3(137);
    z=Altitude;

    % find x and y range of data 
    [Lon, Lat]=ndgrid(ifs.lon, ifs.lat);
    Loc1=[lat1 lon1];

    % get x and y values 
    [x,y] = xy_distance(Lon,Lat,Loc1);
    % find min and max x and y values 
    min_x=min(x(:));  min_y=min(y(:));
    max_x=max(x(:));  max_y=max(y(:));

    % get x and y range for regular distance grid 
    xRange=[min_x-30 max_x+30]; 
    yRange=[min_y-30 max_y+30];
    
    % Make regular distance grid 
    xq=xRange(1):point_spacing(1):xRange(2);
    yq= yRange(1):point_spacing(2): yRange(2);
    [Xq,Yq]=ndgrid(xq,yq);
    
    % Find distance from 0 
    d1=sqrt(Xq.^2+Yq.^2);
    arclen=km2deg(d1);
    
    % Find azimuth
    az=atan2d(Xq,Yq);

    
    % Find lon/lat of the grid points 
    [latout, lonout]=reckon(lat1,lon1,arclen,az);
    
    clear arclen az


    % 3D Gridded Interpolant 
    % reverse order of latitudes 
    [ifs.lat, I]=sort(ifs.lat);
    [lon_grid, lat_grid]=ndgrid(ifs.lon, ifs.lat);
    t=ifs.t(:,I,:,:);

    clear ifs.t ifs.lon ifs.lat
 
    sz=size(t);
    
    % Make 3D grid 
    [Xq2,~,~]=ndgrid(xq,yq,z);
    
    griddedt = nan(size(Xq2));
    
    clear Xq2
    
    for z1 = 1:sz(3)
    
        tslice = t(:,:,z1);
        
        F=griddedInterpolant(lon_grid,lat_grid,tslice, 'linear','none');
    
        griddedt(:,:,z1) = F(lonout,latout);
    
    end
    
    clear t Lon Lat lon_grid lat_grid tslice lonout latout
    
    %griddedt=permute(griddedt,[2,1,3]);
    sz_gt=size(griddedt);

    % find locations of nan values in griddedt
    nan_i_gt=find(isnan(griddedt));

    % replace nan values with median of each horizontal level  
    for m=1:sz_gt(3)
   
        % indicies of nan values in level
        level=squeeze(griddedt(:,:,m));
        nan_m_i=find(isnan(level));
        med=median(level(:),'omitnan');
        level(nan_m_i)=med;
        griddedt(:,:,m)=level;
          
    end
       
    % smooth IFS as AIRS in the horizontal 
    FWHM=[13.5 13.5];
    sigma=(FWHM./[point_spacing(1) point_spacing(2)])./(2*sqrt(2*log(2)));
    griddedt1=imgaussfilt(griddedt,sigma,'padding','replicate','filterdomain', 'spatial');

    % put nan values back in array 
    griddedt1(nan_i_gt)=NaN;


    clear griddedt level nan_m_i

    %% Put IFS onto AIRS grid 

    % find x and y distance for AIRS granule 
    [x_a,y_a] = xy_distance(airs.Data.l1_lon,airs.Data.l1_lat,Loc1);

    sz_airs=size(airs.Data.l1_lon);
    sz_level=size(ifs.level);
     
    temp = nan(sz_airs(1),sz_airs(2),sz_level(1));
    
    
    for z1 = 1:length(ifs.level)
    
        tslice = griddedt1(:,:,z1);
        
        F=griddedInterpolant(Xq,Yq,tslice, 'linear','none');
    
        temp(:,:,z1) = F(x_a,y_a);
    
    end
    
    clear ifs.level tslice Xq Yq



    %% Interpolate onto regular altitude grid
   
    % swap dimensions
    temp=permute(temp,[3,1,2]);
    size(temp);
    
    clear griddedt
    
    % over interpolate
    zq2=HeightRange(1):point_spacing(3):HeightRange(2);
    
    % interpolate 
    grid_t=interp1(z,temp,zq2,'linear');
    size(grid_t);
    clear temp
    
    griddedt1=permute(grid_t,[2,3,1]);
    size(griddedt1);
    clear grid_t

    sz_gt1=size(griddedt1);

  
%%  Replace NaN values and find locations 
   

    % load vertical resolution with altitude 
    airs_v_r=load('airs_vert_res.mat');
    % load AIRS noise with altitude
    airs_noise=load('airs_noise.mat');

    % airs altitudes
    airs_alt=airs_v_r.VertRes.alt;

    % find Airs altitudes within HeightRange 
    alt_i=find(airs_v_r.VertRes.alt>=26 & airs_v_r.VertRes.alt<=55);
    alt=airs_v_r.VertRes.alt(alt_i);
  
    sz_t=size(airs.Data.ret_temp);
    
    % empty array for closest AIRS altitude
    n_airs_alt=[];

    % nan value locations griddedt1
    nan_i_gt1=find(isnan(griddedt1));
    
    % Replace NaNs with median value of each horizontal level 
    % and find closest AIRS altitude to each level
    for k=1:sz_gt1(3)

        % indicies of nan values in level
        level=squeeze(griddedt1(:,:,k));
        nan_m_i=find(isnan(level));
        med=median(level(:),'omitnan');
        level(nan_m_i)=med;
        griddedt1(:,:,k)=level;
        
        % closest airs altitude index to level altitude 
        % minimum difference
        [min_d, min_d_i]=min(abs(airs_alt-zq2(k)));
        % nearest airs altitude
        n_airs_alt=[n_airs_alt airs_alt(min_d_i)];
       
    end


    clear gt1_airs_alt


    %% Smooth to AIRS vertical resolution 
  
    % empty array for altitudes 
    t_ifs=NaN(sz_t(2),sz_t(3),sz_gt1(3));

    TIME=(ifs.ifs_time)*ones(size(airs.Data.l1_lat));
    day_night = which_airs_retrieval(airs.Data.l1_lon, airs.Data.l1_lat,TIME,-1);
    day_frac= nnz(day_night)/numel(day_night);

    midlat_points=(airs.Data.l1_lat >= 30) & (airs.Data.l1_lat < 60);
    midlat_frac= nnz(midlat_points)/numel(midlat_points);

    tropic_points=(airs.Data.l1_lat < 30);
    tropic_frac= nnz(tropic_points)/numel(tropic_points);

    polar_points=(airs.Data.l1_lat >= 60);
    polar_frac= nnz(polar_points)/numel(polar_points);

    % Find resolutions and noise in height range 

    % tropics day
    if  day_frac >= 0.5 && tropic_frac >= 0.5
        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.tropics.day(alt_i); 

        % airs noise folder
        folder='trop_day';

    % tropics night
    elseif day_frac < 0.5 && tropic_frac >= 0.5 

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.tropics.night(alt_i); 

        % airs noise folder
        folder='trop_night';

    % midlatitude day granules
    elseif day_frac >= 0.5 && midlat_frac >= 0.5

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.midlat.day(alt_i); 

        % airs noise folder
        folder='midlat_day';
    
    % midlatitude night
    elseif day_frac < 0.5 && midlat_frac >= 0.5

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.midlat.night(alt_i); 

        % airs noise folder
        folder='midlat_night';

    
    % polar day (winter)
    elseif day_frac >= 0.5 && polar_frac >= 0.5

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.polar_winter.day(alt_i); 

        % airs noise folder
        folder='polar_day';


    % polar night (winter)
    elseif day_frac < 0.5 && polar_frac >= 0.5

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.polar_winter.night(alt_i); 
        % find AIRS noise in height range
        airs_noise=airs_noise.Noise.polar_winter.night(alt_i);

        % airs noise folder
        folder='polar_night';

    end

    % empty array to be filled with very low resolution altitudes to ignore (filled with
    % noise)
    ifs_noise_i=[];

    % loop through resolutions in Height  
    for j=1:length(alt_i)

        % vertical resolution
        v_res=airs_v_res(j);

        % noise
        airs_n=airs_noise(j);

        % find levels in IFS nearest to AIRS altitude
        airs_altitude=alt(j);
        alt_i_ifs=find(n_airs_alt==airs_altitude);

        
        % if v_res does not equal 0
        if v_res~=0
        % calculate sigma from FWHM
        sigma3=(v_res./point_spacing(3))./(2*sqrt(2*log(2)));
        h_sigma=1/2.355;
        %h_sigma=0.01;
        sigma=[h_sigma h_sigma sigma3];

        filt=2*ceil(2*sigma3)+1;

        temp=imgaussfilt3(griddedt1,sigma,'FilterSize',[1 1 filt],'padding','replicate','filterdomain', 'spatial');
         
        % find altitude level 
        temp_level=temp(:,:,alt_i_ifs);

        
        else 
            % if resolution is zero fill levels with level median
            % find median of levels with infinite resolution 
            med_level=median(griddedt1(:,:,alt_i_ifs),[1 2]);
            % make array of medians the same size as the levels
            temp_level=repmat(med_level,sz_gt1(1),sz_gt1(2));

            % save indicies of AIRS altitudes with poor resolution
            ifs_noise_i=[ifs_noise_i j];

        end


        t_ifs(:,:,alt_i_ifs)=temp_level;
        
    end

    clear airs_v_r airs_v_res airs_noise temp

    % put nan values back 
    t_ifs(nan_i_gt1)=NaN;

   

    %% Interpolate to AIRS altitudes 

    t_ifs=permute(t_ifs,[3,1,2]);
    size(t_ifs);
    
    z_t_ifs=zq2;
    zq_airs=alt;

    % interpolate 
    grid_ifs=interp1(z_t_ifs,t_ifs,zq_airs,'linear');
    size(grid_ifs);

    clear t_ifs
    
    t_ifs_a=permute(grid_ifs,[2,3,1]);
    size(t_ifs_a);

    clear grid_ifs

   
    %% Put IFS as AIRS on a regular distance grid 

    % IFS as AIRS structure with NetCDF file names 
    ifs_airs.ret_temp=permute(t_ifs_a,[3,1,2]); 
    ifs_airs.l1_lon=airs.Data.l1_lon;
    ifs_airs.l1_lat=airs.Data.l1_lat;
    ifs_airs.ret_z=alt;
    ifs_airs.l1_time=repmat(ifs.ifs_time,size(airs.Data.l1_lat));

    clear ifs.ifs_time

   
    AirsIn=ifs_airs;


    % use prep AIRS 3D to put on regular distance grid and find Tp 
    [IFS,Spacing,Error, ErrorInfo]  = prep_airs_3d(0,0,'InputStruct',AirsIn,'PreSmooth',[1 1 1],'FullDataDir','./');

    %% Add AIRS noise 
    
    % choose random airs granule with noise out of 30
    file_num=randi(30,1);

    % day/night, region folder
    folder_path=strcat(d_airs_noise,'\',folder);
    FolderInfo=dir(folder_path);
    noise_granules=char(FolderInfo.name);
    noise_granules(1:2,:)=[];

    f_airs=fullfile(folder_path,noise_granules(file_num,:));


    airs_n=load(f_airs);
    % add noise to the IFS
    IFS.Tp = IFS.Tp + airs_n.Tp;
    % add noise to the temperature 
    IFS.T_noise= IFS.ret_temp + airs_n.Tp;

    % smooth IFS.Tp 3x3 boxcar
    IFS.Tp=smoothn(IFS.Tp,[3 3 1]);
    
    %% 2D+1 S-Transform

    sz_a=size(IFS.Tp);
    Point_Spacing=Spacing;

    % nan values in data
    nan_tp=find(isnan(IFS.Tp));

    % if <1/4 of points in granule are NaN fill with median of horizontal
    % level
    if sum(isnan(IFS.Tp(:))) <= (0.25*numel(IFS.Tp))
        for l=1:sz_a(3)
            
            level=IFS.Tp(:,:,l); 

            % find median of level
            med=median(level(:),'omitnan');

            % find locations of NaNs
            nan_i=find(isnan(level));
            % replace nans with mean value 
            level(nan_i)=med;

            IFS.Tp(:,:,l)=level;

        end
    end 

    IN=IFS.Tp; 

    % 2D+1 S-Transform 
    ST=nph_2dst_plus1(IN,nfreqs,Point_Spacing,c,'minwavelengths', minwavelengths, 'maxwavelengths', maxwavelengths);
    
    % if amplitude is not all NaNs
    if ~isnan(max(ST.A(:)))

        ST.z=ifs_airs.ret_z;
        ST.time=ifs.ifs_time;
        ST.spacing=Point_Spacing;
        ST.lon=airs.Data.l1_lon; 
        ST.lat=airs.Data.l1_lat;
    
        A=ST.A; % Wave amplitude 
        Tb = IFS.BG;
        Temperature= IFS.ret_temp;
        P=alt2pres_complex(IFS.ret_z);
        P1=permute(P,[2 3 1]);
        Pressure=repmat(P1,[90 135 1]);
        AirDensity = cjw_airdensity(Pressure,Temperature);
        k = ST.F1;
        l = ST.F2;
        m = ST.F3;
        g = 9.807;
        Nb = 0.02;
        
        ST.MF= (AirDensity./2).*((sqrt((k.^2)+(l.^2)))./abs(m)).*((g/Nb).^2).*((A./Tb).^2);
    
        clear AirDensity P P1 Pressure Tb Temperature k l m
       
        % save nan value locations 
        ST.nan=nan_tp;
        % save index of levels filled with noise 
        ST.ifs_noise_i=ifs_noise_i; 
        % save temperature
        ST.T=IFS.ret_temp;
        ST.T_noise=IFS.T_noise;
        ST.BG=IFS.BG;
    
    
        % remove fields
        fields = {'C','HA','HR','BoostFactor'};
        ST = rmfield(ST,fields);
    
    
        % save S-Transform as structure 
        save(strcat('D:/IFS_AIRS_ST_airs_noise_day_night/IFS_AIRS_ST_',num2str(Day),'_',num2str(h_airs),'h',num2str(m_airs),'m_',num2str(GN),'.mat'),'-struct','ST')
        
    end
    end

end