%% resample ERA5 as AIRS, using ERA5 data from every 3 hours 

% Functions written by others included in this MATLAB script
%%% getnet, nph_2dst_plus1 : Neil Hindley. University of Bath, UK.(nh351@bath.ac.uk)
%%% ecmwf_prs_v3, prep_airs_3d, alt2pres_complex, which_airs_retrieval, 
% cjw_airdensity: Corwin Wright. University of Bath, UK. (cw785@bath.ac.uk) 
%%% smoothn: Anil Gannepali. 25/NOV/2001 (https://uk.mathworks.com/matlabcentral/fileexchange/725-smoothn)

% ERA5 files for AIRS granules
d_e='C:\Users\ejl45\OneDrive - University of Bath\ERA5\ERA5_AIRS_3h';
DirectoryInfo=dir(d_e);
filenames=char(DirectoryInfo.name);
filenames(1:2,:)=[]; 

% AIRS files directory 
a='C:\Users\ejl45\OneDrive - University of Bath\AIRS\2018\';

% AIRS noise directory
d_airs_noise='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise new';

HeightRange=[26 55];

c=[0.25 0.25 0.25];
minwavelengths= [60 60 6];
maxwavelengths= [800 800 45];
nfreqs=500; 
v_point_spacing=0.1;

% loop through files in directory and interpolate onto AIRS grid 
for i=1:length(filenames)

    f_e=fullfile(d_e, strcat(filenames(i,:)));

    era5=load(f_e);

    disp(i);

    Day=era5.airs_day;
    
    % if day is <= 10
    if Day<=10
        disp('>=11')
    Day_str=num2str(Day,'%02d');
    GN=era5.GN; 

    %  find day of year 
    date=datetime(2018, 11, Day);
    d=day(date,'dayofyear');

    % convert matlab time to days and hours
    T_era5=datestr(era5.era5_time,'dd-mm HH:MM:SS');
    h_era5=T_era5(7:8);
    m_era5=T_era5(10:11);
    S_era5=T_era5(13:14);

    % convert airs matlab time to days and hours
    T_airs=datestr(era5.airs_time,'dd-mm HH:MM:SS');
    h_airs=T_airs(7:8);
    m_airs=T_airs(10:11);
    S_airs=T_airs(13:14);


    % AIRS directory 
    d_a=strcat('C:\Users\ejl45\OneDrive - University of Bath\AIRS\2018\', num2str(d));

    GN=num2str(GN);
    % add 0s in front of granule number 
    if numel(GN)==1
       GN=strcat('00',GN);
    end
    if numel(GN)==2
        GN=strcat('0',GN);
    end 

    filename=strcat('airs_2018_',num2str(d),'_', GN,'.nc');
    f_a=fullfile(d_a, filename);
    
    % load corresponding AIRS granule 
    airs=getnet(f_a);

    %% Interpolate ERA5 onto AIRS grid

    % reverse order of latitudes
    [lat, I]=sort(era5.lat);

    % make a grid of lon/lat values
    [lon_grid, lat_grid]=ndgrid(era5.lon, lat);

    t=era5.t(:,I,:);

    sz_airs=size(airs.Data.l1_lon);
    sz_level=size(era5.level);
     
    griddedt = nan(sz_airs(1),sz_airs(2),sz_level(1));
    
    clear Xq2
    
    for z1 = 1:length(era5.level)
    
        tslice = t(:,:,z1);
        
        F=griddedInterpolant(lon_grid,lat_grid,tslice, 'linear','none');
    
        griddedt(:,:,z1) = F(airs.Data.l1_lon,airs.Data.l1_lat);
    
    end
    

    % Altitude of levels 
    [~,z] = ecmwf_prs_v3(137);
    
    %% Interpolate to regular altitude

    % swap dimensions
    temp=permute(griddedt,[3,1,2]);
    size(temp);
    
    clear griddedt
    
    zq2=HeightRange(1):v_point_spacing:HeightRange(2);
    
    % interpolate 
    grid_t=interp1(z,temp,zq2,'linear');
    size(grid_t);
    clear temp
    
    griddedt1=permute(grid_t,[2,3,1]);
    size(griddedt1);
    clear grid_t


    sz_gt1=size(griddedt1);

%% Find closest AIRS altitude to ERA5 model levels

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

    % nan value loactions griddedt1
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
        n_airs_alt=[n_airs_alt airs_alt(min_d_i)];
       
    end

    

    %% Smooth to AIRS vertical resolution 

    % empty array for altitudes 
    t_era5=NaN(sz_t(2),sz_t(3),sz_gt1(3));

    % Find fraction of daytime points and points in midlatitudes, polar
    % regions and tropics
    
    TIME=(era5.era5_time)*ones(size(airs.Data.l1_lat));
    day_night = which_airs_retrieval(airs.Data.l1_lon, airs.Data.l1_lat,TIME,-1);
    day_frac= nnz(day_night)/numel(day_night);

    midlat_points=(airs.Data.l1_lat >= 30) & (airs.Data.l1_lat < 60);
    midlat_frac= nnz(midlat_points)/numel(midlat_points);

    tropic_points=(airs.Data.l1_lat < 30);
    tropic_frac= nnz(tropic_points)/numel(tropic_points);

    polar_points=(airs.Data.l1_lat >= 60);
    polar_frac= nnz(polar_points)/numel(polar_points);

    % tropics day
    if  day_frac >= 0.5 && tropic_frac >= 0.5
        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.tropics.day(alt_i); 
        % find AIRS noise in height range
        airs_noise=airs_noise.Noise.tropics.day(alt_i);

        % airs noise folder
        folder='trop_day';

    % tropics night
    elseif day_frac < 0.5 && tropic_frac >= 0.5 

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.tropics.night(alt_i); 
        % find AIRS noise in height range
        airs_noise=airs_noise.Noise.tropics.night(alt_i);

        % airs noise folder
        folder='trop_night';

    % midlatitude day granules
    elseif day_frac >= 0.5 && midlat_frac >= 0.5

        % find resolutions in height range
        airs_v_res=airs_v_r.VertRes.midlat.day(alt_i); 
        % find AIRS noise in height range
        airs_noise=airs_noise.Noise.midlat.day(alt_i);

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

        % airs noise folder
        folder='polar_night';

    end
    
    % empty array for very low resolution altitudes to ignore (filled with
    % noise)
    era5_noise_i=[];

    % loop through resolutions in Height  
    for j=1:length(alt_i)

        % vertical resolution
        v_res=airs_v_res(j);
        % noise
        airs_n=airs_noise(j);

        % find levels in ERA5 nearest to AIRS altitude 
        airs_altitude=alt(j);
        alt_i_era5=find(n_airs_alt==airs_altitude);
        
        
        % if v_res does not equal 0
        if v_res~=0
        % calculate sigma from FWHM
        sigma3=(v_res./v_point_spacing)./(2.*sqrt(2.*log(2)));
        %sigma3=0.5;

        h_sigma=1/2.355;
        %h_sigma=0.01;
        sigma=[h_sigma h_sigma sigma3];

        filt=2*ceil(2*sigma3)+1;


        temp=imgaussfilt3(griddedt1,sigma,'FilterSize',[1 1 filt],'padding','replicate','filterdomain', 'spatial');


        % add filtered levels closest to airs altitude to array  
        temp_level=temp(:,:,alt_i_era5);
              
              
        else 

            % if resolution value is zero (indicates very low resolutio) 
            % fill levels with level median

            % find median of levels with poor resolution 
            med_level=median(griddedt1(:,:,alt_i_era5),[1 2]);
            % make array of medians the same size as the levels
            temp_level=repmat(med_level,sz_gt1(1),sz_gt1(2));

            % save altitudes with poor resolution
            era5_noise_i=[era5_noise_i j];

           
        end

        % add noise to altitude levels (can be commented out to resample
        % ERA5 as AIRS without adding AIRS noise)
        temp_level=temp_level+(randn(size(temp_level)).*airs_n);


        if any(temp_level<0,'all')
            disp('error 1')
            disp(j)
            disp(i)
        end

        t_era5(:,:,alt_i_era5)=temp_level;

        
    end

    % put nan values back 
    t_era5(nan_i_gt1)=NaN; 

    if any(t_era5<0,'all')
            disp('error 1?')
            disp(j)
            disp(i)
    end

   
    %% Interpolate to AIRS altitudes 

    t_era5=permute(t_era5,[3,1,2]);
    size(t_era5);
    
    z_t_era5=zq2;
    zq_airs=alt;

    % interpolate 
    grid_era5=interp1(z_t_era5,t_era5,zq_airs,'linear');
    size(grid_era5);
    
    t_era5_a=permute(grid_era5,[2,3,1]);
    size(t_era5_a);

     

    %% ERA5 as AIRS structure with NetCDF file names 

    era5_airs.ret_temp=permute(t_era5_a,[3,1,2]); 
    era5_airs.l1_lon=airs.Data.l1_lon;
    era5_airs.l1_lat=airs.Data.l1_lat;
    era5_airs.ret_z=alt;
    era5_airs.l1_time=repmat(era5.era5_time,size(airs.Data.l1_lat));

    Era5In=era5_airs;

    % use prep AIRS 3D to put on regular distance grid and find Tp 
    [Era5,Spacing,Error, ErrorInfo]  = prep_airs_3d(0,0,'InputStruct',Era5In,'PreSmooth',[1 1 1],'FullDataDir','./');



   %% Add AIRS noise 
    
    % choose random AIRS granule containing only noise out of 30 for 
    % region (tropics, midlatitudes, or polar), for daytime and nighttime data 

    file_num=randi(30,1);

    % day/night, region folder
    folder_path=strcat(d_airs_noise,'\',folder);
    FolderInfo=dir(folder_path);
    noise_granules=char(FolderInfo.name);
    noise_granules(1:2,:)=[];

    f_airs=fullfile(folder_path,noise_granules(file_num,:));
 
    airs_n=load(f_airs);
    % add noise to ERA5
    Era5.Tp= Era5.Tp + airs_n.Tp;
    % add noise to the temperature data 
    Era5.T_noise= Era5.ret_temp+ airs_n.Tp; 


    % smooth with 3x3 boxcar in the horizontal 
    Era5.Tp=smoothn(Era5.Tp,[3 3 1]);

   

%%

    sz_a=size(Era5.Tp);
    point_spacing=Spacing;

    % nan values in data
    nan_tp=find(isnan(Era5.Tp));

    % if <1/4 of points in granule are NaN fill with median of horizontal
    % level and fill Airs.ret_temp
    if sum(isnan(Era5.Tp(:))) <= (0.25*numel(Era5.Tp))
        for l=1:sz_a(3)
            
            level=Era5.Tp(:,:,l); 
            level2=Era5.ret_temp(:,:,l);

            % find median of level
            med=median(level(:),'omitnan');
            med2=median(level2(:),'omitnan');

            % find locations of NaNs
            nan_i=find(isnan(level));
            nan_i2=find(isnan(level2));
            % replace nans with median value 
            level(nan_i)=med;
            level2(nan_i2)=med2;

            Era5.Tp(:,:,l)=level;
            Era5.ret_temp(:,:,l)=level2;

        end

    end 

    IN=Era5.Tp; 


    % 2D+1 S-Transform 
    ST=nph_2dst_plus1(IN,nfreqs,point_spacing,c,'minwavelengths', minwavelengths, 'maxwavelengths', maxwavelengths);
    
    % if arrays are not empty
    if ~isnan(max(ST.A(:)))

        ST.z=era5_airs.ret_z;
        ST.time=era5.era5_time;
        ST.lon=airs.Data.l1_lon; 
        ST.lat=airs.Data.l1_lat; 
        
        A=ST.A; % Wave amplitude
        Tb = Era5.BG;
        Temperature= Era5.ret_temp;
        P=alt2pres_complex(Era5.ret_z);
        P1=permute(P,[2 3 1]);
        Pressure=repmat(P1,[90 135 1]);
        AirDensity = cjw_airdensity(Pressure,Temperature);
        k = ST.F1;
        l = ST.F2;
        m = ST.F3;
        g = 9.807;
        Nb = 0.02;
        
        ST.MF= (AirDensity./2).*((sqrt((k.^2)+(l.^2)))./abs(m)).*((g/Nb).^2).*((A./Tb).^2);
        
        if any(AirDensity<0,'all')
        disp('AirDensity error')
        disp(i)
        end
        
        clear Nb g AirDensity P1 Pressure Tb Temperature k l m
        
        % save nan value locations 
        ST.nan=nan_tp;
        % save index of levels filled with noise 
        ST.era5_noise_i=era5_noise_i; 
        ST.T=Era5.ret_temp;
        ST.T_noise=Era5.T_noise; 
        ST.BG=Era5.BG;
        
        
        %save S-Transform as structure 
        save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5\ERA5_AIRS_ST_3h_airs_noise_day_night\ERA5_AIRS_ST_',Day_str,'_',num2str(h_airs),'h',num2str(m_airs),'m_',num2str(GN),'.mat'),'-struct','ST')

    end
    end
end 

