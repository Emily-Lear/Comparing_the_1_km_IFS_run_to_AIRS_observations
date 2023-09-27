%% 2D+1 S-Transform of IFS 

% Functions written by others included in this function: 
%%% nph_2dst_plus1: Neil Hindley. University of Bath, UK. 2017 (nh351@bath.ac.uk)
%%% cjw_airdensity: Corwin Wright. University of Bath, UK.
% (cw785@bath.ac.uk)

Year=2018;
DayRange=1:14;
Month=11;
Directory='/beegfs/scratch/user/u/ejl45/1km_data/'; 
HeightRange=[24 55]; 
xRange=[-4700 4700];
yRange=[-2200 2200];


%% Loop through 1 km data 
% in same locations as AIRS granules and find S-Transform 

c=[0.25 0.25 0.25];
minwavelengths= [60 60 6];
maxwavelengths= [800 800 45];
nfreqs=1000; 
point_spacing=[15 15 1 1];
time_i=1;


% time step
time_step= 180:180:21600;

% loop through timesteps 

for i = 1:length(time_step) 
    
    ts=time_step(i);
    disp(i)


    % convert timestep to string 
    str_ts=num2str(ts);
    
    % add zeros in front of timestep

    % add 2 zeros for 3 digit numbers 
    if numel(str_ts)==3
        str_ts=strcat(['00', str_ts]);
    end
    % add one zero for 4 digit numbers 
    if numel(str_ts)==4
        str_ts=strcat(['0', str_ts]);
    end

    lon1=94;
    lat1=52;

    disp(ts)

    % file name 
    f=strcat(Directory,'ICMSHh3f7+0',str_ts,'_gg_interpGBA_TenthDeg_NH.nc');
    time=ncread(f,'time');

   
    % Convert to serial date number 
    date=double(time)/24 + datenum('1900-01-01 00:00:00');
    Time=datestr(date,'dd-mmm-yyyy HH:MM');
    Hour=Time(13:14);
    Day=Time(1:2);
    
    % get temperature perturbation 
    [ERA5] = prep_era5_3d(f,lon1,lat1, HeightRange, xRange, yRange, time_i, point_spacing);
    
    point_spacing= point_spacing(1:3);
    IN=ERA5.Tp;
    
    % 2D+1 S-Transform 
    ST=nph_2dst_plus1(IN,nfreqs,point_spacing,c,'minwavelengths', minwavelengths, 'maxwavelengths', maxwavelengths);

    % Remove fields from ST 
    r_fields= {'Options', 'Ph', 'C', 'HR'};
    ST=rmfield(ST,r_fields);

   % Calculate momentum flux 

    A=ST.A; %Tp
    Tb= ERA5.Tb; 
    Temperature=ERA5.T;
    Pressure=ERA5.P;
    AirDensity = cjw_airdensity(Pressure,Temperature);
    k = ST.F1;
    l =ST.F2;
    m = ST.F3; 
    g = 9.807;
    Nb = 0.02;
     
    ST.MF= (AirDensity./2).*((sqrt(k.^2+l.^2))./abs(m)).*((g/Nb).^2).*((A./Tb).^2);
    
    clear g Nb Tb Temperature AirDensity Tp k l m Pressure Tp

    % add other fields to ST 
    ST.lon= ERA5.lon;
    ST.lat=ERA5.lat;
    ST.T= ERA5.T;
    ST.P=squeeze(ERA5.P(1,1,:));
    ST.x=ERA5.x;
    ST.y=ERA5.y;
    ST.z=ERA5.z;
    ST.Day=Day; 

    clear ERA5


    % get fieldnames in ST
    fields=fieldnames(ST);
    
    % change precision from double to single 
    for k = [2, 3,6, 8:21] 
    
        field=char(fields(k));
        ST.(field)=single(ST.(field));
    
    end 
    
    clear fields 
    
    % save S-Transform 
    save(strcat('/beegfs/scratch/user/u/ejl45/1km_data_ST_2018/2D+1_','1km_data_',num2str(Day),'_',num2str(Hour),'_',str_ts,'_','ST.mat'),'-struct', 'ST');


end 
