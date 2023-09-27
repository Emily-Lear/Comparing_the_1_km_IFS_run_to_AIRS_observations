% find 10 AIRS noise granules in region and detrend
% for datime and nighttime data

% f='D:/AIRS/2017/305/airs_2017_305_001.nc';
% 
% airs=getnet(f);

%Directory='C:/Users/ejl45/OneDrive - University of Bath/AIRS';
Directory='D:/AIRS';
Year=2020;
DayRange=1:14;
Month=11;


lon_lat_range=load('lon_lat_range.mat');
LatRange=[lon_lat_range.min_lat lon_lat_range.max_lat];
LonRange=[lon_lat_range.min_lon lon_lat_range.max_lon];
LatRange2=[0 90];
LonRange2=[0 180];


waves = airs_find_waves2(Year, Month, DayRange,Directory,LatRange, LonRange, LatRange2, LonRange2);

% % % save variables 
save('find_waves_box_2020.mat', '-struct', 'waves');



%%

%Directory='C:\Users\ejl45\OneDrive - University of Bath\AIRS';
Directory='D:\AIRS';

% load waves in region for year
find_waves=load('find_waves_box2.mat');
%find_waves=load('find_waves_box2.mat');

% colormap
[A]=cbrewer('div', 'RdBu', 30);

% get last 30 granules 
names=fieldnames(find_waves);

% loop through fields
for i=1:length(names)

    name=names{i};
    find_noise.(name)=find_waves.(name)(end-29:end);

end

% Day=zeros(1,30);
% GN=zeros(1,30);

% make strcut of empty arrays for day/night regions 
e_struct=struct('DGN',[],'f_p',[]);

list=struct('trop_day',e_struct,'trop_night',e_struct,'midlat_day',e_struct, ...
    'midlat_night',e_struct,'polar_day',e_struct,'polar_night',e_struct);


count=0;
% loop through noise granules and plot 
for iGranule=length(find_waves.GN):-1:1 %length(find_waves.GN)

    %disp(iGranule)

    % load granule
    i =find_waves.Day(iGranule);
    j =find_waves.GN(iGranule);
    Day_str=num2str(i,'%02d'); 
    GN_str=num2str(j,'%03d');

    % get temperature perturbation (2017)
    [Airs, Spacing] = prep_airs_3d_h(datenum(2018,11,i),j,'PreSmooth', [1 1 1],'KeepOldTime', true, 'HeightRange',[25 55], 'fulldatadir',Directory); 

    % find x and y distances from a lon/lat coordinate 
    lon1=94;
    lat1=52;
    Loc1=[lat1 lon1];
    
    Lon= Airs.l1_lon;
    Lat= Airs.l1_lat;
    
    [x,y] = xy_distance(Lon,Lat,Loc1);
    Airs.x=x;
    Airs.y=y; 

    % get fraction of points in region 
    x_i=find(Airs.x >=-4500 & Airs.x <=4500);
    y_i=find(Airs.y >=-2000 & Airs.y <=2000);
    % find indicies in x and y range
    xy_i=ismember(x_i, y_i);
    % no. points in region
    num_points=sum(xy_i);
    frac_points=num_points/(90*135);

    % if over half the points are in the region
    if frac_points > 0

        % separate granules into day/night and polar/midlat/tropics
        day_night = which_airs_retrieval(Airs.l1_lon, Airs.l1_lat,Airs.l1_time,-1);
        day_night=double(day_night);
        day_frac= nnz(day_night)/numel(day_night);
    
        midlat_points=(Airs.l1_lat >= 30) & (Airs.l1_lat < 60);
        midlat_frac= nnz(midlat_points)/numel(midlat_points);
    
        tropic_points=(Airs.l1_lat < 30);
        tropic_frac= nnz(tropic_points)/numel(tropic_points);
    
        polar_points=(Airs.l1_lat >= 60);
        polar_frac= nnz(polar_points)/numel(polar_points);

        % tropics day
        if  day_frac >= 0.5 && tropic_frac >= 0.5
        
            if frac_points > 0.1
            DayGN=[Day_str,' ',GN_str];
            list.trop_day.DGN=[list.trop_day.DGN string(DayGN)];
            list.trop_day.f_p=[list.trop_day.f_p frac_points];

            
            
%             figure
%             Airs.Tp(:,:,5)=smoothn(Airs.Tp(:,:,5),[3 3]);
%             pcolor(Airs.x, Airs.y, Airs.Tp(:,:,5)); shading flat
%             title(strcat(Day_str, " ", GN_str))
%             colormap(A)
%             caxis([-3 3])
%             axis equal
%             xlim([-5000 5000])
%             ylim([-2500 2500]) 
            end
                
           
        % tropics night
        elseif day_frac < 0.5 && tropic_frac >= 0.5 
            if frac_points > 0.1

            DayGN=[Day_str,' ',GN_str];
            list.trop_night.DGN=[list.trop_night.DGN string(DayGN)];
            list.trop_night.f_p=[list.trop_night.f_p frac_points];

          

            end

        % midlatitude day granules
        elseif day_frac >= 0.5 && midlat_frac >= 0.5

        if frac_points > 0.1
            DayGN=[Day_str,' ',GN_str];
            list.midlat_day.DGN=[list.midlat_day.DGN string(DayGN)];
            list.midlat_day.f_p=[list.midlat_day.f_p frac_points];
        end 

        % midlatitude night
        elseif day_frac < 0.5 && midlat_frac >= 0.5
        
            if frac_points > 0.1
            DayGN=[Day_str,' ',GN_str];
            list.midlat_night.DGN=[list.midlat_night.DGN string(DayGN)];
            list.midlat_night.f_p=[list.midlat_night.f_p frac_points];

            
            end


        % polar day (winter)
        elseif day_frac >= 0.5 && polar_frac >= 0.5
            
            if frac_points > 0.25


            DayGN=[Day_str,' ',GN_str];
            list.polar_day.DGN=[list.polar_day.DGN string(DayGN)];
            list.polar_day.f_p=[list.polar_day.f_p frac_points];

            

            % store variance 


            end

        % polar night (winter)
        elseif day_frac < 0.5 && polar_frac >= 0.5

            day_night = isdaytime(Airs.l1_lat, Airs.l1_lon,Airs.l1_time,'lars');
            day_night=double(day_night);
            day_frac= nnz(day_night)/numel(day_night);
    
            if day_frac >= 0.5 
            disp(day_frac)
            end
            
            if frac_points > 0.25
            DayGN=[Day_str,' ',GN_str];
            list.polar_night.DGN=[list.polar_night.DGN string(DayGN)];
            list.polar_night.f_p=[list.polar_night.f_p frac_points];

            figure
            Airs.Tp(:,:,5)=smoothn(Airs.Tp(:,:,5),[3 3]);
            pcolor(Airs.x, Airs.y, Airs.Tp(:,:,5)); shading flat
            title(strcat(Day_str, " ", GN_str))
            colormap(A)
            caxis([-3 3])
            xlim([-4500 4500])
            ylim([-2000 2000]) 

            figure
            pcolor(Airs.x, Airs.y, double(day_night)); shading flat
            title(strcat(Day_str, " ", GN_str))
            colormap(A)
            caxis([-3 3])
            xlim([-4500 4500])
            ylim([-2000 2000]) 
                    
            
%             xlim([-4500 4500])
%             ylim([-2000 2000]) 

            %pause(4)

            
            end

        end

        

        %pause(4)

    end

end 

% save lists
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_lists.mat'),'-struct', 'list');


%% Remove granules with waves and select 10 with greater 
% point fractions in the region for daytime/nighttime

% file
f='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_lists.mat';
list=load(f);

% new list
list_n=list;
midlat_night_noise=["06 189" "10 003" "10 234" "04 191" "06 238" "07 163" "03 233" ...
    "01 235" "01 169" "02 226" "09 227" "11 225" "14 197" "13 223" "04 244" "03 167" ...
     "08 170" "11 159" "04 174" "14 164" "11 158" "13 173" "11 175" ...
    "14 230" "08 005" "13 190"  "04 158" "01 152" "05 165" ...
    "10 151" "13 157"];
midlat_day_waves=["02 090" "13 087" "03 097" "10 098" "01 050" "06 102" "08 051" "06 053" ...
    "10 049" "13 054" "06 020" "14 028" "03 064" "09 058" "14 061" "10 032" "04 022" ...
    "05 062" "08 018" "13 070" "04 071" "07 027" "08 083" "06 069" "07 076" "01 066" ...
    "08 067" "11 056" "07 060" "09 074" "11 072" "12 063" "10 065" "02 057" "05 095"];
trop_night_waves=["06 222" "08 171" "09 178" "10 169"];
trop_day_waves="13 036";

% day and granule numbers of polar noise
% polar_day_noise=["11 050" "10 043" "08 045" "08 078" "04 082" "12 074" "14 072" "14 083"...
%     "02 063" "13 076" "11 062" "05 068" "09 080" "08 073" "07 066" "12 080" "01 034" ...
%     "08 063" "12 075" "01 047" "02 071" "09 072" "08 065"];
% polar_day_y=[2016 2016 2016 2016 2016 2016 2016 2017 2017 2017 2017 2017 2017 2017 ...
%     2017 2018 2018 2019 2019 2020 2020 2020 2020];
% polar_day_noise=["12 074" "14 072" "09 080" "08 073" "09 064" ...
%      "02 058" "09 072" "11 070" "14 075" "07 074"];
% polar_day_y=[2016 2016 2017 2017 2017 2018 2020 2020 2020 2020];

% polar day new noise (> 0.25 of points in region)
polar_day_noise=[ "14 072" "13 060" "11 062" "09 064" ...
     "07 061" "10 063" "09 072" "10 079" "14 075" "11 070"];
polar_day_y=[2016 2017 2017 2017 2018 2020 ...
    2020 2020 2020 2020];

polar_night_noise=["12 219" "13 226" "11 228" "10 221" "06 225"...
    "07 227" "06 220" "10 229" "09 222" "07 224"];
polar_night_y=[2017 2017 2017 2017 2017 2018 2018 2020 2020 2020];


% remove granules with waves 
ml_day_i=ismember(list_n.midlat_day.DGN,midlat_day_waves);
list_n.midlat_day.DGN(ml_day_i)=[];
list_n.midlat_day.DGN=["05 045" "11 105" "13 103" "06 086" "04 088" ...
    "01 017" "14 060" "11 055" "04 054" "03 047"];

t_night_i=ismember(list_n.trop_night.DGN,trop_night_waves);
list_n.trop_night.DGN(t_night_i)=[];
list_n.trop_night.DGN=["03 217" "10 185" "14 214" "01 219" "09 211" ... 
    "05 215" "10 218" "01 170" "05 182" "11 176"];

t_day_i=ismember(list_n.trop_day.DGN,trop_day_waves);
list_n.trop_day.DGN(t_day_i)=[];
% 10 noise granules 
list_n.trop_day.DGN=["06 035" "10 097" "03 063" "02 039" "12 062"... 
    "13 069" "10 064" "02 105" "06 068" "05 044"];

ml_night_i=~ismember(list_n.midlat_night.DGN,midlat_night_noise);
list_n.midlat_night.DGN(ml_night_i)=[];
list_n.midlat_night.DGN=["13 223" "02 226" "04 191" "06 189" "08 005" ...
    "09 227" "04 158" "10 234" "14 230" "14 197"];

% 1st 15 granules in list
%list_n.trop_day.DGN=list_n.trop_day.DGN(1:15);
%list_n.trop_night.DGN=list_n.trop_night.DGN(1:15);
% 1st 20 granules in list
% list_n.midlat_day.DGN=list_n.midlat_day.DGN(1:30);
% list_n.midlat_night.DGN=list_n.midlat_night.DGN(1:30);

% polar ?
list_n.polar_day.DGN=polar_day_noise;  %.*(0.8/1.85);
list_n.polar_day.Y=polar_day_y;
list_n.polar_night.DGN=polar_night_noise; %*(1.15/1.5);
list_n.polar_night.Y=polar_night_y;

% save new list 
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_list_no_waves.mat'),'-struct', 'list_n');

%% Noise Granules

Directory='C:\Users\ejl45\OneDrive - University of Bath\AIRS';
% colormap
[A]=cbrewer('div', 'RdBu', 30);

% region day/night list directory 
f_list='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_lists.mat';
list=load(f_list);

Day=[6 7 4 5 6 2 10 4 1 8 1 4 9 9 10 10 2 8 4 3 9 7 1 3 5 1 5 9 6 7]';
GN=[7 77 72 79 86 74 82 88 83 100 99 104 107 57 97 48 ...
    56 50 54 47 41 92 33 96 94 49 29 90 52 163]';


% loop through noise granules
for iGranule=1:30

    % load granule
    i = Day(iGranule);
    j = GN(iGranule);
    Day_str=num2str(i,'%02d'); 
    GN_str=num2str(j,'%03d');

    % get temperature perturbation 
    [Airs, Spacing] = prep_airs_3d_h(datenum(2018,11,i),j,'PreSmooth', [1 1 1],'KeepOldTime', true, 'HeightRange',[25 55], 'fulldatadir',Directory); 


    % find x and y distances from a lon/lat coordinate 
            lon1=94;
            lat1=52;
            Loc1=[lat1 lon1];
            
            Lon= Airs.l1_lon;
            Lat= Airs.l1_lat;
            
            [x,y] = xy_distance(Lon,Lat,Loc1);
            Airs.x=x;
            Airs.y=y; 
        
            pcolor(Airs.x, Airs.y, Airs.Tp(:,:,10)); shading flat
            title(strcat(Day_str, " ", GN_str))
            colormap(A)
            caxis([-3 3])
            axis equal
            xlim([-4500 4500])
            ylim([-2000 2000])  
                
            pause(2)
end

DayGN=cat(2, Day, GN);


%% Save AIRS noise granules day / night region

% Day=[6 7 4 5 6 2 10 4 1 8 1 4 9 9 10 10 2 8 4 3 9 7 1 3 5 1 5 9 6 7];
% GN=[7 77 72 79 86 74 82 88 83 100 99 104 107 57 97 48 ...
%     56 50 54 47 41 92 33 96 94 49 29 90 52 163];

% day/granule lists 
f_lists='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_list_no_waves.mat';
lists=load(f_lists);

field_names=["trop_day" "trop_night" "midlat_day" "midlat_night" "polar_day" "polar_night"];

% colormap
[A]=cbrewer('div', 'RdBu', 30);

% Directory
%Directory='C:\Users\ejl45\OneDrive - University of Bath\AIRS';
Directory='D:/AIRS';

%list=lists.midlat_day.DGN;

% loop through regions - day/night
for i=1:length(field_names)
    
    field=field_names(i);

    list=lists.(field).DGN;    

    % loop through granules 
    for j=1:length(list) %5 (day / night switch)
    % 
       D_GN=convertStringsToChars(list(j));
       Day1=str2double(D_GN(1:2));
       GN1=str2double(D_GN(4:6));
       %Day_str=num2str(Day1,'%02d'); 
       %GN_str=num2str(GN1,'%03d');
        
       % year of granule
       if any(1:4==i)
        year=2018;
       elseif i==5 
            year=lists.(field).Y(j);
       elseif i==6
            year=lists.(field).Y(j);
       end

   
       % get temperature perturbation 
       [Airs, Spacing] = prep_airs_3d_h(datenum(year,11,Day1),GN1,'PreSmooth', [1 1 1],'KeepOldTime', true, 'HeightRange',[25 55], 'fulldatadir',Directory,'DayNightFlag',true); 
       
       % if day and night retrieval are used for the granule 
       if ~all(Airs.DayNightFlag(:)) && ~all(Airs.DayNightFlag(:)==0)
           
           day_night=double(Airs.DayNightFlag);
           day_frac= nnz(day_night)/numel(day_night);

           if day_frac >=0.5

               % find rows of nighttime data 

               n_rows=find(all(day_night==0,1));

               if max(n_rows)==135

                  % find indicies of the same number of 
                  % rows of adjacent daytime data 
                  max_d_row=min(n_rows)-1;
                  min_d_row=min(n_rows)-length(n_rows);
                          
               else
                   min_d_row=max(n_rows)+1;
                   max_d_row=max(n_rows)+length(n_rows);
               end

               % flip section of array along track
               d_flipped=flip(Airs.Tp(:,min_d_row:max_d_row,:),2);
               % replace nighttime rows
               Airs.Tp(:,min(n_rows):max(n_rows),:)=d_flipped;

           % mostly nighttime data
           else
               % find row / column indicies of daytime data
               d_rows=find(all(day_night==0,1));
               % find indicies of rows of adjacent nighttime data 
               if max(d_rows)==135
                  max_n_row=min(d_rows)-1;
                  min_n_row=min(d_rows)-length(d_rows);
               else
                   min_n_row=max(d_rows)+1;
                   max_n_row=max(d_rows)+length(d_rows);
               end
               % flip section of array along track
               n_flipped=flip(Airs.Tp(:,min_n_row:max_n_row,:),2);
               % replace daytime rows 
               Airs.Tp(:,min(d_rows):max(d_rows),:)=n_flipped;
           end


       end 

        % find x and y distances from a lon/lat coordinate 
                lon1=94;
                lat1=52;
                Loc1=[lat1 lon1];
                
                Lon= Airs.l1_lon;
                Lat= Airs.l1_lat;
                
                [x,y] = xy_distance(Lon,Lat,Loc1);
                Airs.x=x;
                Airs.y=y; 

                %figure
                
                Airs.Tp_s(:,:,5)=smoothn(Airs.Tp(:,:,5),[3 3]);
                pcolor(Airs.x, Airs.y, Airs.Tp_s(:,:,5)); shading flat
                title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
                colormap(A)
                caxis([-3 3])
                axis equal
%                 xlim([-4500 4500])
%                 ylim([-2000 2000])  
                
                % flip noise granule (across track direction)
                Airs.Tp_2=flip(Airs.Tp,1);
                % along track 
                Airs.Tp_3=flip(Airs.Tp,2);


%                 figure
% 
%                 Airs.Tp_2_s(:,:,5)=smoothn(Airs.Tp_2(:,:,5),[3 3]);
%                 pcolor(Airs.x, Airs.y, Airs.Tp_2_s(:,:,5)); shading flat
%                 title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
%                 colormap(A)
%                 caxis([-1 1])
%                 axis equal
                    
                %pause(3)
                Airs_1.Tp=Airs.Tp;
                Airs_2.Tp=Airs.Tp_2;
                Airs_3.Tp=Airs.Tp_3;
    
               % save Airs.Tp 
               save(strjoin(['C:\Users\ejl45\OneDrive - University of Bath\AIRS\Airs_noise\Airs noise new\',field, '\AIRS_noise_Tp_',num2str(Day1,'%02.f'),'_',num2str(year),'_',num2str(GN1,'%03.f'), '.mat'],''),'-struct','Airs_1');
               % save Airs flipped (across track direction)
               save(strjoin(['C:\Users\ejl45\OneDrive - University of Bath\AIRS\Airs_noise\Airs noise new\',field, '\AIRS_noise_Tp_',num2str(Day1,'%02.f'),'_',num2str(year),'_',num2str(GN1,'%03.f'), '_flipped_XT.mat'],''),'-struct','Airs_2');
               %along track
               save(strjoin(['C:\Users\ejl45\OneDrive - University of Bath\AIRS\Airs_noise\Airs noise new\',field, '\AIRS_noise_Tp_',num2str(Day1,'%02.f'),'_',num2str(year),'_',num2str(GN1,'%03.f'), '_flipped_AT.mat'],''),'-struct','Airs_3');

    end 
end


%% Loop through Polar night and mirror night data 

Directory='D:/AIRS';

% day/granule lists 
f_lists='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\airs noise\day_night_regions_list_no_waves.mat';
lists=load(f_lists);

%list=lists.polar_day.DGN;   
names=["trop_day" "trop_night" "midlat_day" "midlat_night" "polar_day" "polar_night"];

% colormap
[A]=cbrewer('div', 'RdBu', 30);
for i=6 %:length(names)

    name=names(i);
    list=lists.(name).DGN;

    % loop through granules 
    for j=1:length(list)

        D_GN=convertStringsToChars(list(j));
        Day1=str2double(D_GN(1:2));
        GN1=str2double(D_GN(4:6));
        year=lists.(name).Y(j);
        %year=2018;

        % get temperature perturbation 
       [Airs, Spacing] = prep_airs_3d_h(datenum(year,11,Day1),GN1,'PreSmooth', [1 1 1],'KeepOldTime', true, 'HeightRange',[25 55], 'fulldatadir',Directory,'DayNightFlag',true); 
    
        % find x and y distances from a lon/lat coordinate 
                lon1=94;
                lat1=52;
                Loc1=[lat1 lon1];
                
                Lon= Airs.l1_lon;
                Lat= Airs.l1_lat;
                
                [x,y] = xy_distance(Lon,Lat,Loc1);
                Airs.x=x;
                Airs.y=y; 

                % find day/night
                day_night = which_airs_retrieval(Airs.l1_lon, Airs.l1_lat,Airs.l1_time,-1);


                figure
                
                Airs.Tp_s=smoothn(Airs.Tp(:,:,5),[3 3]);
                pcolor(Airs.Tp_s); shading flat % Airs.x,Airs.y,
                title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
                colormap(A)
                caxis([-3 3])
                axis equal

                figure

                day_night=double(day_night);
                pcolor(day_night); shading flat
                title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
                colormap("parula")
                caxis([0 1])
                axis equal
%                 
%                 figure
% %                 % region 
%                 pcolor(Airs.l1_lat)
%                 %caxis([60, 61]); 
%                 shading flat
%                 axis equal
%                 colorbar

%                 figure
%                 day_night_2=double(day_night_2);
%                 pcolor(day_night_2); shading flat
%                 title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
%                 colormap("parula")
%                 caxis([0 1])
%                 axis equal

% Airs.x, Airs.y,
%                 day_night=Airs.DayNightFlag;
%                 pcolor(day_night); shading flat
%                 title(strcat(D_GN(1:2), " ", D_GN(4:6)," ",num2str(year)))
%                 colormap("parula")
%                 caxis([0 1])
%                 axis equal

                pause(2)

    end

end
