%% Plot color maps of AIRS, IFS and ERA5 for each half day, 1-10th Nov 2018
% 3x3 presmoothing, no noise


% Directories 
d_a='C:\Users\ejl45\OneDrive - University of Bath\AIRS\2D+1 ST 3';
d_e='C:\Users\ejl45\OneDrive - University of Bath\ERA5\ERA5_AIRS_ST_3h_airs_noise_day_night';
d_i='D:\IFS_AIRS_ST_airs_noise_day_night';
d_i_s3d='C:\Users\ejl45\OneDrive - University of Bath\1km Data\IFS_as_AIRS_sampled_s3d_airs_noise_day_night';



DI_e=dir(d_e);
filenames_e=char(DI_e.name);
filenames_e(1:2,:)=[];

DI_a=dir(d_a);
filenames_a=char(DI_a.name);
filenames_a(1:2,:)=[];

DI_i=dir(d_i);
filenames_ifs=char(DI_i.name);
filenames_ifs(1:2,:)=[];
% find Day and GN
GN_ifs=filenames_ifs(:,23:25);
Day_ifs=filenames_ifs(:,13:14);
% DayGN
DayGN=cat(2, Day_ifs, GN_ifs);


DI_i_s3d=dir(d_i_s3d);
filenames_s3d=char(DI_i_s3d.name);
filenames_s3d(1:2,:)=[];

% Find indicies of IFS files in AIRS, ERA5 and IFS_s3d

% find indicies of matching day and granule numbers for sampled T
Day_s3d=filenames_s3d(:,21:22);
GN_s3d=filenames_s3d(:,31:33);
% concatenate strings
DayGN_s3d=cat(2, Day_s3d, GN_s3d);
% find indices
[~,s3d_i]=ismember(str2num(DayGN),str2num(DayGN_s3d));

% Find indicies for ERA5
Day_e=filenames_e(:,14:15);
GN_e=filenames_e(:,24:26);
% concatenate strings 
DayGN_e=cat(2, Day_e,GN_e);
[~, e_i]=ismember(str2num(DayGN),str2num(DayGN_e));

% Find indicies for AIRS
Day_a=filenames_a(:,14:15);
GN_a=filenames_a(:,24:26);
% concatenate strings 
DayGN_a=cat(2, Day_a, GN_a);
[~,a_i]=ismember(str2num(DayGN),str2num(DayGN_a));


%% Loop through granules 


data_name_struct = struct('ifs',0,'s3d',1,'era5',2,'airs',3);
DataNames = fieldnames(data_name_struct);

for i=1:length(filenames_ifs) %9 %339

    disp(i)

    % load data granule

    % ifs 
    file_ifs=filenames_ifs(i,:);
    filepath_ifs=fullfile(d_i,file_ifs);
    ifs=load(filepath_ifs);

    % s3d
    s3d_index=s3d_i(i);
    file_s3d=filenames_s3d(s3d_index,:);
    filepath_s3d=fullfile(d_i_s3d,file_s3d);
    s3d=load(filepath_s3d);

    % AIRS 
    a_index=a_i(i);
    file_a=filenames_a(a_index,:);
    filepath_a=fullfile(d_a, file_a);
    airs=load(filepath_a);

    % ERA5 
    e_index=e_i(i);
    file_e=filenames_e(e_index,:);
    filepath_e=fullfile(d_e, file_e);
    era5=load(filepath_e);

    % AIRS granule time
    day=file_a(14:15);
    day_num=str2double(day);

    %if day_num>=11
        disp('day>=11')
    time=file_a(17:22);

    GN=file_a(24:26);

    data=struct('ifs',ifs,'s3d',s3d,'era5',era5,'airs',airs);
    data_names=fieldnames(data);


    % Find x and y distances for granule and resampled models 
    lon1=94;
    lat1=52;
    Loc1=[lat1 lon1];
    
    Lon= double(airs.lon);
    Lat= double(airs.lat);
    
    [x,y] = xy_distance(Lon,Lat,Loc1);
    data.x=x;
    data.y=y; 

    % find x and y distance for ifs 2 
    Lon2= s3d.lon;
    Lat2= s3d.lat;
    
    [x2,y2] = xy_distance(Lon2,Lat2,Loc1);
    data.s3d.x=x2;
    data.s3d.y=y2; 

    % loop through datasets

    for data_i=1:length(data_names)

        name=data_names{data_i};

        % find level at 39 km alt
        alt= 39;
        [~,level_i]=min(abs(data.(name).z-alt));

        LAT=double(data.(name).lat);
        LON=double(data.(name).lon);
        TIME=double(data.(name).time);

        % day / night
        data.(name).DN = which_airs_retrieval(LON,LAT,TIME,-1);

%         % save day/night and region number
%         if strcmp(name,'s3d') || strcmp(name,'era5')
%             data.(name).DN_region=data.(name).DN_region*ones(size(LAT));
% 
%         end


        % wave properties
    
        data.(name).A=squeeze(data.(name).A(:,:,level_i));

        % temperature
        data.(name).T=squeeze(data.(name).T(:,:,level_i));
        
        % Horizontal Wavelength at 39 km altitude 
        data.(name).HW= 1./sqrt(squeeze(data.(name).F1(:,:,level_i)).^2+squeeze(data.(name).F2(:,:, level_i)).^2);
    
        % Vertical Wavelength at 39 km altitude 
        data.(name).VW = 1./squeeze(data.(name).F3(:,:,level_i));

        %% COMPUTE ANGLES AND PROJECT:

        % if all F3 are +ve, change all frequencies to negative 
        p=find(data.(name).F3>0);
        data.(name).F3(p)=-data.(name).F3(p);
        data.(name).F2(p)=-data.(name).F2(p);
        data.(name).F1(p)=-data.(name).F1(p);

        data.(name).kh = hypot(data.(name).F1(:,:,level_i),data.(name).F2(:,:,level_i));
        sz = size(squeeze(data.(name).A(:,:,1)));
        
        % find angle CLOCKWISE from along track direction...
        ang_at = atan2d(data.(name).F1(:,:,level_i),data.(name).F2(:,:,level_i));
        
        % find azimuth of AT direction CLOCKWISE from north...
        xt_mid = floor(sz(1)./2);
        [~,az_at] = distance(data.(name).lat(xt_mid,1:end-1,1),data.(name).lon(xt_mid,1:end-1,1), ...
            data.(name).lat(xt_mid,2:end,1),data.(name).lon(xt_mid,2:end,1));
        az_at(end+1) = az_at(end);
        %az_at = repmat(az_at,size(ST.kh,1),1,size(ST.kh,3));
        az_at = repmat(az_at,size(data.(name).kh,1),1);
        
        % add angles... (since they both should be clockwise from north)
        data.(name).AN = wrapTo180(az_at + ang_at);

        % Angle 
        %data.(name).AN=atan2d(data.(name).F2(:,:,level_i), data.(name).F1(:,:,level_i));
    
        % horizontal momentum flux 
        data.(name).MF=squeeze(data.(name).MF(:,:,level_i));

        % zonal momentum flux
        data.(name).MF_z=squeeze(data.(name).MF_z(:,:,level_i));

        % meridional momentum flux 
        data.(name).MF_m=squeeze(data.(name).MF_m(:,:,level_i));

        if any(data.(name).MF<0,'all')
            disp('error')
        end
        
        data.(name).Tp=data.(name).IN(:,:,level_i);

        
        % find indicies of vertical wavelengths < 6 km & > 45 km in
        % granule and set to NaN
        wave_i=find(data.(name).VW<6 | data.(name).VW>45);
    
        % set points where VW > 45 km to NaN
        data.(name).VW(wave_i)=NaN;
        data.(name).MF(wave_i)=NaN;
        data.(name).MF_z(wave_i)=NaN;
        data.(name).MF_m(wave_i)=NaN;
        data.(name).AN(wave_i)=NaN;
        data.(name).HW(wave_i)=NaN;
        data.(name).A(wave_i)=NaN;
        % no data removed from Tp
        %data.(name).Tp(wave_i)=NaN; 


    
    end

    field_names=["A","HW","VW","MF","AN","Tp","T","MF_z","MF_m"];

    for field_i=1:length(field_names)
    
        field=field_names(field_i);
    
        % double AIRS data 
        data.airs.(field)=double(data.airs.(field));
    end
    
    % save structure 
    save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise day_night\data_ST_',day,'_',time,'_',GN,'.mat'),'-struct', 'data');
        
end
    




%% Find min / max values for variables in half days for plots 

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data 1';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 
% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',4:21, 'd3',22:38, 'd4',39:59, 'd5',60:76, ...
    'd6',77:94, 'd7',95:111, 'd8',112:131, 'd9',132:148, 'd10',149:168, ...
    'd11',169:185, 'd12',186:205, 'd13',206:223, 'd14',224:243, 'd15',244:260, ...
    'd16',261:279, 'd17',280:296, 'd18',297:317, 'd19',318:334, 'd20', ...
    335:352, 'd21',353:370);

half_days=fieldnames(file_indices);

wave_properties=["A" "HW" "VW" "MF" "AN" "Tp"];

% empty arrays for min / max values for each half day 
for fieldi=1:length(wave_properties)
 
    field=wave_properties(fieldi);

    c.airs.min.(field)=[]; c.airs.max.(field)=[];
    c.era5.min.(field)=[]; c.era5.max.(field)=[];
    c.ifs.min.(field)=[];  c.ifs.max.(field)=[];
    c.s3d.min.(field)=[];  c.s3d.max.(field)=[];

end 

% loop through half days
for i=1:length(half_days)

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    for fieldi=1:length(wave_properties)
        
        field=wave_properties(fieldi);

        % empty arrays for max/ min values of granules 
        cb.airs.min.(field)=[]; cb.airs.max.(field)=[];
        cb.era5.min.(field)=[]; cb.era5.max.(field)=[];
        cb.ifs.min.(field)=[];  cb.ifs.max.(field)=[];
        cb.s3d.min.(field)=[];  cb.s3d.max.(field)=[];

    end

    % loop through granules in half day 
    for j=1:length(file_i)
   
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % find points in x and y range 
        xi=find(data.x>-4500 & data.x<4500);
        % find y values in x range 
        y_values=data.y(xi);
        yi=find(y_values>-2000 & y_values<2000);

        % loop through fields
        for fieldi=1:length(wave_properties)
        
            field=wave_properties(fieldi);

            % find average of points in distance range 
            airs.(field)=data.airs.(field)(xi);
            airs.(field)=airs.(field)(yi);
            % find min/max and append to list 
            cb.airs.min.(field)=[cb.airs.min.(field) min(airs.(field)(:))];
            cb.airs.max.(field)=[cb.airs.max.(field) max(airs.(field)(:))];
    
            era5.(field)=data.era5.(field)(xi);
            era5.(field)=era5.(field)(yi);
            % find min/max and append to list
            cb.era5.min.(field)=[cb.era5.min.(field) min(era5.(field)(:))];
            cb.era5.max.(field)=[cb.era5.max.(field) max(era5.(field)(:))];
    
            ifs.(field)=data.ifs.(field)(xi);
            ifs.(field)=ifs.(field)(yi);
            % find min/max and append to list
            cb.ifs.min.(field)=[cb.ifs.min.(field) min(ifs.(field)(:))];
            cb.ifs.max.(field)=[cb.ifs.max.(field) max(ifs.(field)(:))];
    
            s3d.(field)=data.s3d.(field)(xi);
            s3d.(field)=s3d.(field)(yi);
            % find min/max and append to list
            cb.s3d.min.(field)=[cb.s3d.min.(field) min(s3d.(field)(:))];
            cb.s3d.max.(field)=[cb.s3d.max.(field) max(s3d.(field)(:))];

        end
        

    end

    % loop through fields
    for fieldi=1:length(wave_properties)
    
        field=wave_properties(fieldi);

        % find min/ max values for half day and append to lists
        c.airs.min.(field)=[c.airs.min.(field) min(cb.airs.min.(field))];
        c.era5.min.(field)=[c.era5.min.(field) min(cb.era5.min.(field))];
        c.ifs.min.(field)=[c.ifs.min.(field) min(cb.ifs.min.(field))];
        c.s3d.min.(field)=[c.s3d.min.(field) min(cb.s3d.min.(field))];
    
        c.airs.max.(field)=[c.airs.max.(field) max(cb.airs.max.(field))];
        c.era5.max.(field)=[c.era5.max.(field) max(cb.era5.max.(field))];
        c.ifs.max.(field)=[c.ifs.max.(field) max(cb.ifs.max.(field))];
        c.s3d.max.(field)=[c.s3d.max.(field) max(cb.s3d.max.(field))];

    end

    % last AIRS granule in half day 
    filename=filenames_d(index,:);
    day=filename(9:10);
    time=filename(12:17);
    %c.day=str2double(day);
    %c.time=str2double(time);
    

end

% save structure 
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties\min max\min_max.mat'),'-struct', 'c');



   


%% Amplitude 

% remove points below 40th percentile 

% % Percentiles
% f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
% p=load(f_prc);

% day_prc=2:2:20;
% night_prc=3:2:21;

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data no noise new';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 
% min/max file
% f_min_max='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties\min max\min_max.mat';
% cb=load(f_min_max);

% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);

half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data

for i=2:length(half_days) %11

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='compact';
    t.Padding='tight';
    set(gcf,'color','w');


    % loop through granules in half day 
    for j=1:length(file_i)
   
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

%         % smooth granule with 7x7 filter 
%         s.airs.A=smoothn(data.airs.A,[7 7]);
%         s.ifs.A=smoothn(data.ifs.A,[7 7]);
%         s.s3d.A=smoothn(data.s3d.A,[7 7]);
%         s.era5.A=smoothn(data.era5.A,[7 7]);
% 
%         if ismember(i,day_prc)
%             data.airs.A(s.airs.A<p.airs.day)=NaN;
%             data.ifs.A(s.ifs.A<p.ifs.day)=NaN;
%             data.s3d.A(s.s3d.A<p.s3d.day)=NaN;
%             data.era5.A(s.era5.A<p.era5.day)=NaN;
% 
%         elseif ismember(i,night_prc)
%             data.airs.A(s.airs.A<p.airs.night)=NaN;
%             data.ifs.A(s.ifs.A<p.ifs.night)=NaN;
%             data.s3d.A(s.s3d.A<p.s3d.night)=NaN;
%             data.era5.A(s.era5.A<p.era5.night)=NaN;
%         end

%         % remove values below percentiles 
%         data.airs.A(abs(data.airs.A)<0.7890)=NaN;
%         data.ifs.A(abs(data.ifs.A)<0.5815)=NaN;
%         data.s3d.A(abs(data.s3d.A)<0.5764)=NaN;
%         data.era5.A(abs(data.era5.A)<0.5455)=NaN;


        % plot granules  

        % AIRS Amplitude
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, data.airs.A); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'OrRd', 30);
        colormap(ax(1),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.airs.min.A(i) cb.airs.max.A(i)])
        set(gca,'TickDir','out');
        clim([0 3])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4700 4700])
        ylim([-2200 2200])        
        hold(ax(1),'on') 
        
        
        % ERA5 Amplitude
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, data.era5.A); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'OrRd', 30);
        colormap(ax(2),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.era5.min.A(i) cb.era5.max.A(i)])
        clim([0 3])
        set(gca,'TickDir','out');
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 
      
        
        % IFS Amplitude
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, data.ifs.A); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'OrRd', 30);
        colormap(ax(3),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.ifs.min.A(i) cb.ifs.max.A(i)])
        set(gca,'TickDir','out');
        clim([0 3])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(3),'on') 
        
     
        % IFS s3d Amplitude
        %set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.s3d.x, data.s3d.y, data.s3d.A); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'OrRd', 30);
        colormap(CT);
        colormap(ax(4),CT);
        c=colorbar;
        c.Label.String='Amplitude (K)';
        %caxis([cb.s3d.min.A(i) cb.s3d.max.A(i)])
        set(gca,'TickDir','out');
        clim([0 3])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000]) 
        hold(ax(4),'on') 
        

    end

   


    % add topography to plots
    ax(1)=nexttile(1);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(2)=nexttile(2);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(3)=nexttile(3);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(4)=nexttile(4);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 
    

    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['A ',day_1,' - ',day,' ',time_1(1:2),':',time_1(4:5),' - ', ...
         time(1:2),':', time(4:5)]);
    
    
%     % save plot 
%     plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Amplitude\colormap_',day,'_',time,'_',GN,'.png');
%     saveas(gcf,plotname);

end 

%% Horizontal Wavelength

% Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise new 2\';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% % min/max file
% f_min_max='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties\min max\min_max.mat';
% cb=load(f_min_max);

% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);


half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data


for i=2 %:length(half_days) %11

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

        % smooth granule with 7x7 filter 
        s.airs.A=smoothn(data.airs.A,[7 7]);
        s.ifs.A=smoothn(data.ifs.A,[7 7]);
        s.s3d.A=smoothn(data.s3d.A,[7 7]);
        s.era5.A=smoothn(data.era5.A,[7 7]);

        if ismember(i,day_prc)
            data.airs.HW(s.airs.A<p.airs.day)=NaN;
            data.ifs.HW(s.ifs.A<p.ifs.day)=NaN;
            data.s3d.HW(s.s3d.A<p.s3d.day)=NaN;
            data.era5.HW(s.era5.A<p.era5.day)=NaN;

        elseif ismember(i,night_prc)
            data.airs.HW(s.airs.A<p.airs.night)=NaN;
            data.ifs.HW(s.ifs.A<p.ifs.night)=NaN;
            data.s3d.HW(s.s3d.A<p.s3d.night)=NaN;
            data.era5.HW(s.era5.A<p.era5.night)=NaN;
        end

        % AIRS HW
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, data.airs.HW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Blues', 30);
        colormap(ax(1),CT);
        %c=colorbar;
        %c.Label.String='Horizontal Wavelength (km)';
        %caxis([cb.airs.min.HW(i) cb.airs.max.HW(i)])
        set(gca,'TickDir','out');
        clim([60 490])
        %caxis([60 550])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(1),'on') 

        min_airs=min(data.airs.HW(:));

        
        % ERA5 HW
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, data.era5.HW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Blues', 30);
        colormap(ax(2),CT);
        %c=colorbar;
        %c.Label.String='Horizontal Wavelength (km)';
        %caxis([cb.era5.min.HW(i) cb.era5.max.HW(i)])
        set(gca,'TickDir','out');
        %caxis([0 550])
        clim([60 490])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 

        max_era5=max(data.era5.HW(:));


        % IFS HW
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, data.ifs.HW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Blues', 30);
        colormap(ax(3),CT);
        %c=colorbar;
        %c.Label.String='Horizontal Wavelength (km)';
        %caxis([cb.ifs.min.HW(i) cb.ifs.max.HW(i)])
        set(gca,'TickDir','out');
        %caxis([0 550]) %60
        clim([60 490])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        %hold on
        hold(ax(3),'on') 

        min_ifs=min(data.ifs.HW(:));
        

        % IFS s3d HW
        %set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, data.s3d.HW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Blues', 30);
        colormap(CT);
        colormap(ax(4),CT);
        c=colorbar;
        c.Label.String='Horizontal Wavelength (km)';
        %caxis([cb.s3d.min.HW(i) cb.s3d.max.HW(i)])
        set(gca,'TickDir','out');
        %caxis([0 550])
        clim([60 490])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        hold(ax(4),'on') 

        %max_s3d=max(data.s3d.HW(:))
        

    end

    % add topography to plots
    ax(1)=nexttile(1);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(2)=nexttile(2);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(3)=nexttile(3);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(4)=nexttile(4);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 
    
    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['HW ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);


%     % save plot 
%     plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Horizontal Wavelength\colormap_',day,'_',time,'_',GN,'.png');
%     saveas(gcf,plotname);

end 

%% Vertical Wavelength 6

% Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise new 2';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 
% % min/max file
% f_min_max='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties\min max\min_max.mat';
% cb=load(f_min_max);

% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);


half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data


for i=2:length(half_days)

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

        % smooth granule with 7x7 filter 
        s.airs.A=smoothn(data.airs.A,[7 7]);
        s.ifs.A=smoothn(data.ifs.A,[7 7]);
        s.s3d.A=smoothn(data.s3d.A,[7 7]);
        s.era5.A=smoothn(data.era5.A,[7 7]);

        if ismember(i,day_prc)
            data.airs.VW(s.airs.A<p.airs.day)=NaN;
            data.ifs.VW(s.ifs.A<p.ifs.day)=NaN;
            data.s3d.VW(s.s3d.A<p.s3d.day)=NaN;
            data.era5.VW(s.era5.A<p.era5.day)=NaN;

        elseif ismember(i,night_prc)
            data.airs.VW(s.airs.A<p.airs.night)=NaN;
            data.ifs.VW(s.ifs.A<p.ifs.night)=NaN;
            data.s3d.VW(s.s3d.A<p.s3d.night)=NaN;
            data.era5.VW(s.era5.A<p.era5.night)=NaN;
        end

        % plot granules
        
        % AIRS VW
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, data.airs.VW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Greens', 30);
        colormap(ax(1),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.airs.min.VW(i) cb.airs.max.VW(i)])
        set(gca,'TickDir','out');
        caxis([6 45])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(1),'on') 

        
        % ERA5 VW
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, data.era5.VW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Greens', 30);
        colormap(ax(2),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.era5.min.VW(i) cb.era5.max.VW(i)])
        set(gca,'TickDir','out');
        caxis([6 45])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 
        %hold on
        

        % IFS VW
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, data.ifs.VW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Greens', 30);
        colormap(ax(3),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([cb.ifs.min.VW(i) cb.ifs.max.VW(i)])
        set(gca,'TickDir','out');
        caxis([6 45])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        %hold on
        hold(ax(3),'on') 
      

        % IFS s3d VW
        %set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, data.s3d.VW); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Greens', 30);
        colormap(CT);
        colormap(ax(4),CT);
        c=colorbar;
        c.Label.String='Vertical Wavelength (km)';
        %caxis([cb.s3d.min.VW(i) cb.s3d.max.VW(i)])
        set(gca,'TickDir','out');
        caxis([6 45])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        hold(ax(4),'on') 

        min_s3d=min(data.s3d.VW(:));
       

    end

    % add topography to plots
    ax(1)=nexttile(1);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(2)=nexttile(2);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(3)=nexttile(3);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(4)=nexttile(4);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 
    
    set(c,'TickDir','out');
    c.Layout.Tile='east';
%     
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['VW ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);
    
    
    % save plot 
    plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Vertical Wavelength 6\colormap_',day,'_',time,'_',GN,'.png');
    saveas(gcf,plotname);

end 

%% Momentum Flux

% Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise new 2';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 
% min/max file
f_min_max='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties\min max\min_max.mat';
cb=load(f_min_max);

% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);


half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data

for i=1:length(half_days)

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

        % smooth granule with 7x7 filter 
        s.airs.A=smoothn(data.airs.A,[7 7]);
        s.ifs.A=smoothn(data.ifs.A,[7 7]);
        s.s3d.A=smoothn(data.s3d.A,[7 7]);
        s.era5.A=smoothn(data.era5.A,[7 7]);

        if ismember(i,day_prc)
            data.airs.MF(s.airs.A<p.airs.day)=NaN;
            data.ifs.MF(s.ifs.A<p.ifs.day)=NaN;
            data.s3d.MF(s.s3d.A<p.s3d.day)=NaN;
            data.era5.MF(s.era5.A<p.era5.day)=NaN;

        elseif ismember(i,night_prc)
            data.airs.MF(s.airs.A<p.airs.night)=NaN;
            data.ifs.MF(s.ifs.A<p.ifs.night)=NaN;
            data.s3d.MF(s.s3d.A<p.s3d.night)=NaN;
            data.era5.MF(s.era5.A<p.era5.night)=NaN;
        end

        % plot granules
        
        % AIRS MF
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, log10(data.airs.MF)); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Purples', 30);
        colormap(ax(1),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([log10(cb.airs.min.MF(i)) log10(cb.airs.max.MF(i))])
        set(gca,'TickDir','out');
        caxis([-5 -1.5])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(1),'on') 

        
        % ERA5 MF
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, log10(data.era5.MF)); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Purples', 30);
        colormap(ax(2),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([log10(cb.era5.min.MF(i)) log10(cb.era5.max.MF(i))])
        set(gca,'TickDir','out');
        caxis([-5 -1.5])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 
        

        % IFS MF
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, log10(data.ifs.MF)); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Purples', 30);
        colormap(ax(3),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([log10(cb.ifs.min.MF(i)) log10(cb.ifs.max.MF(i))])
        set(gca,'TickDir','out');
        caxis([-5 -1.5])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(3),'on') 
        

        % IFS s3d MF
        %set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, log10(data.s3d.MF)); shading flat 
        % change colormap 
        CT=cbrewer('seq', 'Purples', 30);
        colormap(CT);
        colormap(ax(4),CT);
        c=colorbar;
        c.Label.String='log10(Momentum Flux)';
        %caxis([log10(cb.s3d.min.MF(i)) log10(cb.s3d.max.MF(i))])
        set(gca,'TickDir','out');
        caxis([-5 -1.5])
        %set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        min_s3d=min(log10(data.s3d.MF(:)));

        hold(ax(4),'on') 
        

    end

    % add topography to plots
    ax(1)=nexttile(1);
%     c=colorbar;
%     c.Label.String='Vertical Wavelength (km)'; 
%     set(c,'TickDir','out')
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(2)=nexttile(2);
%     c=colorbar;
%     c.Label.String='Vertical Wavelength (km)'; 
%     set(c,'TickDir','out')
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(3)=nexttile(3);
%     c=colorbar;
%     c.Label.String='Vertical Wavelength (km)'; 
%     set(c,'TickDir','out')
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(4)=nexttile(4);
%     c=colorbar;
%     c.Label.String='Vertical Wavelength (km)'; 
%     set(c,'TickDir','out')
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 
    
    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['MF ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);
    
    % save plot 
    plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Momentum Flux\colormap_',day,'_',time,'_',GN,'.png');
    saveas(gcf,plotname);

end 


%%  Angle

% Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise new 2';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 
% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);


% % make stucture of indexes for each half day
% file_indices=struct('d1',1:3, 'd2',4:21, 'd3',22:38, 'd4',39:59, 'd5',60:76, ...
%     'd6',77:94, 'd7',95:111, 'd8',112:131, 'd9',132:148, 'd10',149:168, ...
%     'd11',169:185, 'd12',186:205, 'd13',206:223, 'd14',224:243, 'd15',244:260, ...
%     'd16',261:279, 'd17',280:296, 'd18',297:317, 'd19',318:334, 'd20', ...
%     335:352, 'd21',353:370);

half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data


for i=2:length(half_days)

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

        % smooth granule with 7x7 filter 
        s.airs.A=smoothn(data.airs.A,[7 7]);
        s.ifs.A=smoothn(data.ifs.A,[7 7]);
        s.s3d.A=smoothn(data.s3d.A,[7 7]);
        s.era5.A=smoothn(data.era5.A,[7 7]);

        if ismember(i,day_prc)
            data.airs.AN(s.airs.A<p.airs.day)=NaN;
            data.ifs.AN(s.ifs.A<p.ifs.day)=NaN;
            data.s3d.AN(s.s3d.A<p.s3d.day)=NaN;
            data.era5.AN(s.era5.A<p.era5.day)=NaN;

        elseif ismember(i,night_prc)
            data.airs.AN(s.airs.A<p.airs.night)=NaN;
            data.ifs.AN(s.ifs.A<p.ifs.night)=NaN;
            data.s3d.AN(s.s3d.A<p.s3d.night)=NaN;
            data.era5.AN(s.era5.A<p.era5.night)=NaN;
        end

        % plot granules
        
        % AIRS AN
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, data.airs.AN); shading flat 
        % change colormap 
        [A]=cbrewer('div', 'BrBG', 30);
        colormap(A);
        %colormap default;
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([-180 180])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(1),'on') 

        
        % ERA5 AN
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, data.era5.AN); shading flat 
        % change colormap 
        [A]=cbrewer('div', 'BrBG', 30);
        colormap(A);
        %colormap default;
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([-180 180])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 
        

        % IFS AN
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, data.ifs.AN); shading flat 
        % change colormap 
        [A]=cbrewer('div', 'BrBG', 30);
        colormap(A);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([-180 180])
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(3),'on') 
       

        % IFS s3d AN
        %set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, data.s3d.AN); shading flat 
        % change colormap 
        [A]=cbrewer('div', 'BrBG', 30);
        colormap(A);
        c=colorbar;
        c.Label.String=['Angle (',char(176),')'];
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        caxis([-180 180])
        set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        hold(ax(4),'on') 
   

    end

%     %add topography to plots
%     ax(1)=nexttile(1);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(2)=nexttile(2);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(3)=nexttile(3);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(4)=nexttile(4);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 
    
    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['AN ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);
    
    % save plot 
    plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Angle\colormap_',day,'_',time,'_',GN,'.png');
    saveas(gcf,plotname);

end 

%% Tp - airs noise
% remove amplitudes below percentiles 

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise day_night';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% % Percentiles
% f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc_50.mat';
% p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;
 
% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);

% 163 ->
half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data

for i=2:length(half_days) %11

    figure
    set(gcf,'WindowState','maximized')

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)

        %disp(j)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end
        
%         % smooth granule with 7x7 filter 
%         s.airs.A=smoothn(data.airs.A,[7 7]);
%         s.ifs.A=smoothn(data.ifs.A,[7 7]);
%         s.s3d.A=smoothn(data.s3d.A,[7 7]);
%         s.era5.A=smoothn(data.era5.A,[7 7]);
% 
%         if ismember(i,day_prc)
%             data.airs.Tp(abs(s.airs.A)<p.airs.day)=NaN;
%             data.ifs.Tp(abs(s.ifs.A)<p.ifs.day)=NaN;
%             data.s3d.Tp(abs(s.s3d.A)<p.s3d.day)=NaN;
%             data.era5.Tp(abs(s.era5.A)<p.era5.day)=NaN;
% 
%         elseif ismember(i,night_prc)
%             data.airs.Tp(abs(s.airs.A)<p.airs.night)=NaN;
%             data.ifs.Tp(abs(s.ifs.A)<p.ifs.night)=NaN;
%             data.s3d.Tp(abs(s.s3d.A)<p.s3d.night)=NaN;
%             data.era5.Tp(abs(s.era5.A)<p.era5.night)=NaN;
%         end

%         % remove values below percentiles 
%         data.airs.Tp(abs(data.airs.Tp)<0.7890)=NaN;
%         data.ifs.Tp(abs(data.ifs.Tp)<0.5815)=NaN;
%         data.s3d.Tp(abs(data.s3d.Tp)<0.5764)=NaN;
%         data.era5.Tp(abs(data.era5.Tp)<0.5455)=NaN;

        % plot granules
        
        % AIRS Tp
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, data.airs.Tp); shading flat 
        % change colormap 
        CT=cbrewer('div', 'RdBu', 30);
        colormap(ax(1),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        %caxis([-2.5 2.5])
        caxis([-3 3])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
%         xlim([-4700 4700])
%         ylim([-2200 2200])
        hold(ax(1),'on') 

        
        % ERA5 Tp
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, data.era5.Tp); shading flat 
        % change colormap 
        CT=cbrewer('div', 'RdBu', 30);
        colormap(ax(2),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        %caxis([-2.5 2.5])
        caxis([-3 3])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5 1')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
%         xlim([-4700 4700])
%         ylim([-2200 2200])
        hold(ax(2),'on') 
        

        % IFS Tp
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, data.ifs.Tp); shading flat 
        % change colormap 
        CT=cbrewer('div', 'RdBu', 30);
        colormap(ax(3),CT);
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        %caxis([-2.5 2.5])
        caxis([-3 3])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
%         xlim([-4700 4700])
%         ylim([-2200 2200])
        %hold on
        hold(ax(3),'on') 
        

        % IFS s3d Tp
        set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, data.s3d.Tp); shading flat 
        % change colormap 
        CT=cbrewer('div', 'RdBu', 30);
        colormap(CT);
        colormap(ax(4),CT);
        c=colorbar;
        c.Label.String='Temperature Perturbation (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        %caxis([-2.5 2.5])
        caxis([-3 3])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
%         xlim([-4700 4700])
%         ylim([-2200 2200])
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        hold(ax(4),'on') 
       

    end

    % add topography to plots
    ax(1)=nexttile(1);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(2)=nexttile(2);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(3)=nexttile(3);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 

    ax(4)=nexttile(4);
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
    hold on
    contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 

    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['Tp ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);
    
    
    % save plot 
%     plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\ERA5_IFS_AIRS_colormaps_new\Tp\colormap_',day,'_',time,'_',GN,'_HR.png');
%     saveas(gcf,plotname);

end 


%% Day / Night

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data no noise new';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];
 

% make stucture of indexes for each half day
file_indices=struct('d1',1:3, 'd2',5:21, 'd3',22:42, 'd4',44:61, 'd5',64:82, ...
    'd6',85:102, 'd7',104:123, 'd8',125:144, 'd9',145:162, 'd10',165:182, ...
    'd11',185:203, 'd12',205:224, 'd13',225:244, 'd14',248:265, 'd15',268:286, ...
    'd16',288:306, 'd17',307:327, 'd18',329:346, 'd19',350:368, 'd20', ...
    371:387, 'd21',389:404);


half_days=fieldnames(file_indices);

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[30 30 1 1];
% add topography 
[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% loop through half days of data

for i=2:length(half_days)

    figure

    % get file indexes 
    field=half_days{i};
    file_i=file_indices.(field);

    t=tiledlayout(2,2);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    for j=1:length(file_i)
        
        % file index
        index=file_i(j);
        % load file 
        f=fullfile(D,filenames_d(index,:));
        data=load(f);

        % save time for 1st granule
        if j==1

         filename_1=filenames_d(index,:);
         day_1=filename_1(9:10);
         time_1=filename_1(12:17);   

        end

        % plot granules
        
        % AIRS Day / Night
        %set(gcf,'color','w');
        ax(1)=nexttile(1);
        pcolor(data.x, data.y, double(data.airs.DN)); shading flat 
        % change colormap 
        colormap(ax(1),'parula(2)');
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([0 1])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('AIRS')
        %xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(1),'on') 
        
        %T_era5=datestr(data.era5.time,'dd-mm HH:MM:SS')
        
        % ERA5 Day / Night
        %set(gcf,'color','w');
        ax(2)=nexttile(2);
        pcolor(data.x, data.y, double(data.era5.DN_region)); shading flat 
        % change colormap 
        colormap(ax(2),'parula(6)');
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([1 6])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('ERA5')
        %xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        hold(ax(2),'on') 

        % IFS Day / Night
        %set(gcf,'color','w');
        ax(3)=nexttile(3);
        pcolor(data.x, data.y, double(data.ifs.DN)); shading flat 
        % change colormap 
        colormap(ax(3),'parula(2)');
        %c=colorbar;
        %c.Label.String='Amplitude (K)';
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([0 1])
        %c.Ticks=0:1;
        %c.TickLabels=["Night" "Day"];
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 1')
        xlabel('Distance (km)')
        ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        %hold on
        hold(ax(3),'on') 
        

        % IFS s3d Day / Night
        set(gcf,'color','w');
        ax(4)=nexttile(4);
        pcolor(data.x, data.y, double(data.s3d.DN_region)); shading flat 
        % change colormap 
        colormap(ax(4),'parula(6)');
        c=colorbar;
        %c.Label.String='Day / Night';
%         c.Ticks=0:1;
%         c.TickLabels=["Night" "Day"];
        %caxis([min(e_ST.A(:)) max(e_ST.A(:))])
        set(gca,'TickDir','out');
        caxis([1 6])
        %set(gca,'TickDir','out');
        % set(c,'TickDir','out');
        % equal axes 
        axis equal
        title('IFS 2')
        xlabel('Distance (km)')
        %ylabel('Distance (km)')
        xlim([-4500 4500])
        ylim([-2000 2000])
        
        hold(ax(4),'on') 

%         % print times 
%         j
%         T_s3d=datestr(data.s3d.time,'dd-mm HH:MM:SS')
    
       

    end

%     % add topography to plots
%     ax(1)=nexttile(1);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(2)=nexttile(2);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(3)=nexttile(3);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 
% 
%     ax(4)=nexttile(4);
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0], 'LineWidth', 1.5);
%     hold on
%     contour(Topo_km.x, Topo_km.y, Topo_km.elev', 0:1:500,'-','color',[0.5 0.5 0.5]); 


    hold off 

    set(c,'TickDir','out');
    c.Layout.Tile='east';
    
     % AIRS granule time
     filename=filenames_d(index,:);
     day=filename(9:10);
     time=filename(12:17);
    
     GN=filename(19:21);
    
     title(t, ['D / N ',day_1,' ',time_1(1:2),':',time_1(4:5),' - ', ...
        day,' ', time(1:2),':', time(4:5)]);
    
    
%     % save plot 
%     plotname=strcat('C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps_new\Day Night\colormap_',day,'_',time,'_',GN,'.fig');
%     saveas(gcf,plotname,'fig');

end 

%% Plot Day/night for lon lat graphs 


lat=0:0.1:90;
lon=0:0.1:180;
[LON, LAT]=ndgrid(lon, lat);
TIME=datenum(2018,11,1,6,0,0);

% find x and y distances 
lon1=94;
lat1=52;
Loc1=[lat1 lon1];

[x,y] = xy_distance(LON,LAT,Loc1);

% day / night
OUT = which_airs_retrieval(L0N,LAT,TIME,-1);

figure
pcolor(x, y, double(OUT)); shading flat
axis equal
xlim([-4500 4500])
ylim([-2000 2000])

figure 
pcolor(x,y, LON); shading flat
axis equal
clim([0 180])
title('Lon')
colorbar
xlim([-4500 4500])
ylim([-2000 2000])

figure
pcolor(x,y,LAT); shading flat
axis equal
clim([0 90])
title('Lat')
colorbar
xlim([-4500 4500])
ylim([-2000 2000])

