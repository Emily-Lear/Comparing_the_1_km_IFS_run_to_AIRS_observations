%% Find 2D+1 S - Transform for AIRS

% Functions written by others included in this MATLAB script:
%%% prep_airs_3d_h (same function as prep_airs_3d: Corwin Wright. University of Bath, UK. (cw785@bath.ac.uk), but with altitude range
% as an optional input)
%%% alt2pres_complex, cjw_airdensity: Corwin Wright. University of Bath, UK. (cw785@bath.ac.uk) 
%%% smoothn: Anil Gannepali. 25/NOV/2001 (https://uk.mathworks.com/matlabcentral/fileexchange/725-smoothn)
%%% nph_2dst_plus1: Neil Hindley. University of Bath, UK.(nh351@bath.ac.uk)

Year=2018;
DayRange=1:14;
Month=11;
Directory='C:\Users\ejl45\OneDrive - University of Bath\AIRS';


lon_lat_range=load('lon_lat_range.mat');
LatRange=[lon_lat_range.min_lat lon_lat_range.max_lat];
LonRange=[lon_lat_range.min_lon lon_lat_range.max_lon];
LatRange2=[0 90];
LonRange2=[0 180];
 

%% Uncomment this section to save AIRS granules in region, sorted in order of variance

%[waves] = airs_find_waves(Year, Month, DayRange,Directory,LatRange, LonRange, LatRange2, LonRange2);

%save variables 
%save('find_waves_box2.mat', '-struct', 'waves');

%% Loop through all Airs granules to get S-Transform 

c=[0.25 0.25 0.25];
minwavelengths= [60 60 6];
maxwavelengths= [800 800 45];
nfreqs=500;

find_waves=load('find_waves_box2.mat');

for iGranule = 1:length(find_waves.Day)
    
    i = find_waves.Day(iGranule);
    j = find_waves.GN(iGranule);
    Day_str=num2str(i,'%02d'); 
    GN_str=num2str(j,'%03d');

 
    disp(j);
        
    % get temperature perturbation 
    [Airs, Spacing] = prep_airs_3d_h(datenum(2018,11,i),j,'PreSmooth', [1 1 1],'KeepOldTime', true, 'HeightRange',[25 55], 'fulldatadir',Directory); 


    % smooth using a 3x3 boxcar 
    Airs.Tp=smoothn(Airs.Tp,[3 3 1]);
    
    time=mean(Airs.OriginalTime(:));
    T = datetime(time,'ConvertFrom','epochtime','Epoch','2000-01-01');
    T=datestr(T,'dd-mm HH:MM:SS');
    Hour=T(7:8);
    Minute=T(10:11);
    
    point_spacing=Spacing;
    IN=Airs.Tp;

    % 2D+1 S-Transform 
    ST=nph_2dst_plus1(IN,nfreqs,point_spacing,c,'minwavelengths', minwavelengths, 'maxwavelengths', maxwavelengths);


    % Remove fields from ST 
    r_fields= {'Options', 'Ph', 'C', 'HR'};
    ST=rmfield(ST,r_fields);

    % Calculate momentum flux 

    A=ST.A; % Wave amplitude
    Tb = Airs.BG;
    Temperature= Airs.ret_temp;
    P=alt2pres_complex(Airs.ret_z);
    P1=permute(P,[2 3 1]);
    Pressure=repmat(P1,[90 135 1]);
    AirDensity = cjw_airdensity(Pressure,Temperature);
    k = ST.F1;
    l = ST.F2;
    m = ST.F3;
    g = 9.807;
    Nb = 0.02;
    
    ST.MF= (AirDensity./2).*((sqrt((k.^2)+(l.^2)))./abs(m)).*((g/Nb).^2).*((A./Tb).^2);
    

    clear Nb g AirDensity P1 Pressure Tb Temperature k l m

    ST.T= Airs.ret_temp;
    ST.Tp=Airs.Tp;
    ST.P=P;
    ST.lon= Airs.l1_lon;
    ST.lat= Airs.l1_lat;
    ST.time=Airs.l1_time;
    ST.z=Airs.ret_z;
    ST.Spacing=Spacing;
    ST.Day=5; 
    ST.OriginalTime=Airs.OriginalTime;
    ST.BG=Airs.BG;


    % get fieldnames in ST
    fields=fieldnames(ST);

    % change precision from double to single 
    for k = [2, 3,6, 8:22] 
    
        field=char(fields(k));
        ST.(field)=single(ST.(field));

    end 

    clear fields 

    % save S-Transform 
    save(strcat('C:/Users/ejl45/OneDrive - University of Bath/AIRS/2D+1 ST 3/2D+1_ST_AIRS_',Day_str,'_',num2str(Hour),'h',num2str(Minute),'m_',GN_str,'.mat'),'-struct', 'ST');
    
end 

