%% Get S-Transform for ERA5 and save

% Functions written by others included in this function: 
%%% nph_2dst_plus1: Neil Hindley. University of Bath, UK. 2017 (nh351@bath.ac.uk)
%%% cjw_airdensity: Corwin Wright. University of Bath, UK.
% (cw785@bath.ac.uk)

Year=2018;
Month=11;
Directory= 'D:/ERA5_2018';
HeightRange=[26 55];
xRange=[-4700 4700];
yRange=[-2200 2200];

c=[0.25 0.25 0.25];
minwavelengths= [60 60 6];
maxwavelengths= [800 800 45];
%nfreqs=1000; 
nfreqs=500; %300
point_spacing=[30 30 1 1];

find_waves=load('find_waves.mat');
[x, y]=meshgrid(0:2:180, 0:2:90);
store.N=zeros(size(x));
store.A=store.N;

for Day = 1:10

    lon1=94;
    lat1=52;

    disp(Day)

    % find day of year 
    t=datetime(Year, Month, Day);
    DayofYear=day(t,'dayofyear');
  
    f=strcat(Directory,'/era5_',num2str(Year),'d',num2str(DayofYear),'.nc');
    time=ncread(f,'time');

    for time_i=1:length(time)
       
        time1=time(time_i);

        % Convert to serial date number 
        date=double(time1)/24 + datenum('1900-01-01 00:00:00');
        Time=datestr(date,'dd-mmm-yyyy HH:MM');
        Hour=Time(13:14);
    
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

        A=ST.A;
        Tb= ERA5.Tb; 
        Temperature=ERA5.T;
        Pressure=ERA5.P;
        AirDensity = cjw_airdensity(Pressure,Temperature);
        k = ST.F1;
        l =ST.F2;
        m = ST.F3; 
        g = 9.807;
        Nb = 0.02;
        
        ST.MF= (AirDensity./2).*((sqrt((k.^2)+(l.^2)))./abs(m)).*((g/Nb).^2).*((A./Tb).^2);
        
        clear g Nb Tb Temperature AirDensity Tp k l m 

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
        for k = [2, 3, 6, 8:21] 
        
            field=char(fields(k));
            ST.(field)=single(ST.(field));
        
        end 
        
        clear fields 


        % save S-Transform 
        save(strcat('C:/Users/ejl45/OneDrive - University of Bath/ERA5/2D+1 ST/2D+1_','ERA5_',num2str(ST.Day),'_',num2str(Hour),'_','ST.mat'),'-struct', 'ST');
    end 

end 

