%% Example maps 

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% % make stucture of indexes for each half day
file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
    'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
    'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
    'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
    116:121, 'd21',123:127);

half_days=fieldnames(file_indices);


data_names=["airs" "era5" "ifs" "s3d"];
names=["AIRS" "ERA5" "IFS 1" "IFS 2"];

fields=["Tp" "A" "HW" "VW" "MF" "AN"];
%field_names=["Tp" "Amplitude" "Horizontal Wavelength"]
cb_labels=["Tempertature perturbation (K)" "Amplitude (K)" "Horizontal Wavelength (km)" ...
    "Vertical wavelength (km)" "log10( Momentum flux )" ['Angle (',char(176),')']];

% Colormaps
CT_1=cbrewer('div', 'RdBu', 30);
CT_2=cbrewer('seq', 'OrRd', 30);
CT_3=cbrewer('seq', 'Blues', 30);
CT_4=cbrewer('seq', 'Greens', 30);
CT_5=cbrewer('seq', 'Purples', 30); 
CT_6=cbrewer('div', 'BrBG', 30);
    
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5,'CT_6',CT_6);
CT_names = fieldnames(CT_s);

cb_l_limits=[-3 0 60 12 -5 -180];
cb_u_limits=[3 3 490 35 -1.5 180];

%granules=[3];


% get file indexes 
half_day=half_days{3}; %11 %13 %17 %19
file_i=file_indices.(half_day);
    
% loop through wave properties 
for field_i=1:length(fields)
    
    field=fields(field_i);
    CT_name=CT_names{field_i};
    CT=CT_s.(CT_name);
    cb_ll=cb_l_limits(field_i);
    cb_ul=cb_u_limits(field_i);
    
    f=figure; 
    %set(gcf,'WindowState','maximized')
    
    t=tiledlayout(1,4);
    t.TileSpacing='tight';
    t.Padding='tight';
    set(gcf,'color','w');

    

    % loop through datasets 
    for data_i=1:length(data_names)
    
        name=data_names(data_i);

        % loop through granules in half day 
        for i=4 %:length(granules) %1:length(file_i)
        
            % file index
            index=file_i(i);
            %index=granules(i);
            % load file 
            file=fullfile(D,filenames_d(index,:));
            data=load(file);

            % log10( Momentum Flux )
            if strcmp(field,'MF')
                data.(name).MF=log10(data.(name).MF);
            end
        
            % save time for 1st granule
            if isfield(data,'gn')
    
                day_1=char(data.gn.day(1));
                time_1=char(data.gn.airs_time(1)); 
                gn_1=char(data.gn.GN(1));
            end

            
            ax(data_i)=nexttile(data_i);
            pcolor(data.x, data.y, data.(name).(field)); shading flat 
            % change colormap 
            colormap(ax(data_i),CT);
            set(gca,'TickDir','out');
            clim([cb_ll cb_ul])
            % equal axes 
            axis equal
            title(names(data_i))
            if data_i==1 %|| data_i==3
               ylabel('Distance (km)')
            end
            if data_i==3 || data_i==4 || data_i==1 || data_i==2
               xlabel('Distance (km)')
            end            
            xlim([-4700 4700])
            %xlim([-1500 2000]) % d11 %d19
            %xlim([-1000 3000]) %d13
            xlim([-2500 1500]) % d17 %d3
            ylim([-2000 2000])             
            set(gca,'Color',[0.8 0.8 0.8])
            set(gca,'TickDir','out');
            hold(ax(data_i),'on') 

        end
        
    end
    c=colorbar;
    c.Label.String=cb_labels(field_i);
    set(c,'TickDir','out');
    c.Layout.Tile='east';

    f.Position(1)=f.Position(1)-250;
    f.Position(2)=f.Position(2)-40;
    f.Position(3)=f.Position(3)*2;
    f.Position(4)=f.Position(4)*0.7;
    
    % AIRS granule time
    day=char(data.gn.day(end));
    time=char(data.gn.airs_time(end));
    GN=char(data.gn.GN(end));
    
    
    title(t, strjoin([field,' ',time_1(1:2),':',time_1(4:5),' ',day_1,' - ', ...
     time(1:2),':', time(4:5),' ',day],''));

    pause(1)
    
    
%   % save plot 
    exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps' ...
        '\example maps\d3\colormap_',field,'_',day,'_',time,'.png'],''),'Resolution',300);
  
    
end

%% Example maps % 09 (A, HW, VW, MF)

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% % make stucture of indexes for each half day
file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
    'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
    'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
    'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
    116:121, 'd21',123:127);

half_days=fieldnames(file_indices);


data_names=["airs" "ifs" "s3d" "era5"];
names=["AIRS" "IFS 1" "IFS 2" "ERA5 1"];

fields=["A" "HW" "VW" "MF_z" "MF_m"];
%field_names=["Tp" "Amplitude" "Horizontal Wavelength"]
cb_labels=["Amplitude (K)" "Horizontal Wavelength (km)" ...
    "Vertical wavelength (km)" "Zonal MF (mPa)" "Meridional MF (mPa)"]; %"log10( Momentum flux (Pa) )"
field_names=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" ...
    "Zonal MF" "Meridional MF"]; %"Momentum Flux"

labels=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)" "(m)" "(n)" "(o)" "(p)" ...
    "(q)" "(r)" "(s)" "(t)"];
labels=reshape(labels,4,5);

% Colormaps
%CT_1=cbrewer('div', 'RdBu', 30);
CT_1=cbrewer('seq', 'OrRd', 30);
CT_2=cbrewer('seq', 'Blues', 30);
CT_3=cbrewer('seq', 'Greens', 30);
CT_4_a=flipud(cbrewer('div', 'BrBG',24)); % MF_z
CT_5_a=cbrewer('div', 'PRGn',24); % MF_m
CT_4=ones(30,3);
CT_4(1:12,:)=CT_4_a(1:12,:);
CT_4(19:30,:)=CT_4_a(13:24,:);
CT_5=ones(30,3);
CT_5(1:12,:)=CT_5_a(1:12,:);
CT_5(19:30,:)=CT_5_a(13:24,:);
    
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5);
CT_names = fieldnames(CT_s);

cb_l_limits=[0 60 12 -4 -4]; % MF: -5, MF_z & m: -5
cb_u_limits=[3 490 35 4 4]; % MF: -1.5, MF_z & m: 5

%granules=[3];

tile=1:20; %16
tile=reshape(tile,4,5);

% get file indexes 
half_day=half_days{11}; %11 %19
file_i=file_indices.(half_day);

f=figure;

f.Position(1)=f.Position(1)-250;
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.7; %*0.4; %width %0.84
f.Position(4)=1170*0.8; %0.63

t=tiledlayout(5,4);
t.TileSpacing='tight';
%t.Padding='tight';
set(gcf,'color','w');
    
% loop through wave properties 
for field_i=1:length(fields) % Tp - MF
    
    field=fields(field_i);
    CT_name=CT_names{field_i};
    CT=CT_s.(CT_name);
    cb_ll=cb_l_limits(field_i);
    cb_ul=cb_u_limits(field_i);

    t_col=tile(:,field_i);
    l_col=labels(:,field_i);
    
        

    % loop through datasets 
    for data_i=1:length(data_names)
    
        name=data_names(data_i);

        % loop through granules in half day 
        for i=3 %:length(granules) %1:length(file_i)
        
            % file index
            index=file_i(i);
            %index=granules(i);
            % load file 
            file=fullfile(D,filenames_d(index,:));
            data=load(file);

            % log10( Momentum Flux )
            % if strcmp(field,'MF')
            %     data.(name).MF=log10(data.(name).MF);

            
            % end
            if strcmp(field,'MF_z') || strcmp(field,'MF_m')
                
                data.(name).(field)=-data.(name).(field)*10^3;
                
                % log MF_z and MF_m
                MF=data.(name).(field);
                MF_sign=sign(MF);

                % log absolute values + 0.005
                log_MF=log10(abs(MF)+0.005);
                min_log_MF=min(log_MF(:));
                % make all values positive
                log_MF=log_MF+abs(log10(0.005)); 

                %log_MF(zero_i)=0;
                % change signs back 
                data.(name).(field)=(log_MF.*MF_sign);

                % colorbar values
                cb_labels_2=[-1e2 -10 -1 -0.1 0 0.1 1 10 1e2];
                cb_ticks=[];

                for i=1:9

                    cbl=cb_labels_2(i);
                    cb_ticks(i)=(log10(abs(cbl)+0.005)+abs(log10(0.005)))*sign(cbl);

                end

                %     %max_MF_z=max(data.(name).MF_z(:))
                % if strcmp(field,'MF_z')
                %     min_MF_z=min(data.(name).MF_z(:))
                %     max_MF_z=max(data.(name).MF_z(:))
                % elseif strcmp(field,'MF_m')
                %     min_MF_m=min(data.(name).MF_m(:))
                %     max_MF_m=max(data.(name).MF_m(:))
                % end

            end
        
            % save time for 1st granule
            if isfield(data,'gn')

    
                day_1=char(data.gn.day(1));
                time_1=char(data.gn.airs_time(1)); 
                gn_1=char(data.gn.GN(1));
            end
            t_col(data_i);
            
            ax(t_col(data_i))=nexttile(t_col(data_i));
            pcolor(data.x, data.y, data.(name).(field)); shading flat 
            % change colormap 
            colormap(ax(t_col(data_i)),CT);
            set(gca,'TickDir','out');
            clim([cb_ll cb_ul])
            
            % equal axes 
            axis equal
            if field_i==1
                t2=title(names(data_i));
            end   
            if data_i==1
               
                st=subtitle(field_names(field_i));

                st.Rotation=90;
                st.Position=[-3100 0]; %-2900
                st.FontWeight='bold';

               % row_t=text(-3200,0,field_names(field_i),'HorizontalAlignment','center');
               % row_t.Rotation=90;
               % row_t.FontWeight='bold';
               % %row_t.Position=[-2000 -2000];

               ylabel('Distance (km)')               
            else
                set(gca,'Yticklabel',[])
            end
            
            if field_i==5
               xlabel('Distance (km)')
            else
                set(gca,'Xticklabel',[])
            end  
            label=l_col(data_i);
            t_l=text(0,0,label);
            t_l.Units='normalized';
            t_l.Position=[0.02 0.93];
            t_l.FontSize=8.5;
            ax(t_col(data_i)).FontSize=7.5;
            xlim([-4700 4700])
            xlim([-1500 2000]) % d11 %d19
            %xlim([-1000 3000]) %d13
            %xlim([-2500 1500]) % d17 %d3
            ylim([-2000 2000])             
            set(gca,'Color',[0.8 0.8 0.8])
            set(gca,'TickDir','out');
            if data_i==4
                c=colorbar;
                
                c.Layout.Tile='east';
                if field_i==4 || field_i==5 
                    c.Ticks=cb_ticks;
                    c.TickLabels=cb_labels_2;
                end
             
            end
            
        end
        
    end

    c.Position(1)=0;
    c.Position(3)=c.Position(3)*0.5;
    c.Label.String=cb_labels(field_i);
    set(c,'TickDir','out');

    
   
end


% AIRS granule time
day=char(data.gn.day(end));
time=char(data.gn.airs_time(end));
GN=char(data.gn.GN(end));

t1=title(t, [time_1(1:2),':',time_1(4:5),' ',' - ', ...
    time(1:2),':', time(4:5),' ','UTC ',day(2),'^{th}',' November']);
t1.FontSize=9;

%fontsize(9,"points");
    
% save plot 
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps' ...
    '\example maps\tiled layout\colormap_',day,'_',time,'_log_MF.png'],'Resolution',300);


%% Example maps Tp (airs noise and no noise added)

% Directory - with AIRS noise
d='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night';
DI_d=dir(d);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% directory - no noise added
d_n='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes no noise day_night\';
DI_n=dir(d_n);
filenames_n=char(DI_n.name);
filenames_n(1:2,:)=[];

% topography filepath
f_t='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\topography\topo';
Topo_km=load(f_t);

% map1=[linspace(0,1,25)' linspace(0.35,1,25)' linspace(0.1,1,25)'];
% map2=[linspace(1,0.22,25)' linspace(1,0.14,25)' linspace(1,0.01,25)'];
map1=[linspace(0,0,10)' linspace(0.25,0.35,10)' linspace(0.005,0.01,10)'];
map2=[linspace(0,1,120)' linspace(0.35,1,120)' linspace(0.01,1,120)'];
map3=[linspace(1,0.05,220)' linspace(1,0.05,220)' linspace(1,0.05,220)'];

% map2(1,:)=[];
% map=[map1; map2];
map2(1,:)=[];
map3(1,:)=[];
map=[map1; map2; map3];

% % make stucture of indexes for each half day
file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
    'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
    'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
    'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
    116:121, 'd21',123:127);

half_days=fieldnames(file_indices);
% get file indexes 
half_day=half_days{11}; %11 %19
file_i=file_indices.(half_day);

data_names=["airs" "ifs" "s3d" "era5"]; % "airs"
names=["AIRS" "IFS 1" "IFS 2" "ERA5 1"]; % "AIRS"

fields=["A" "HW" "VW" "MF"];
%field_names=["Tp" "Amplitude" "Horizontal Wavelength"]
labels=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)"];
labels=reshape(labels,4,2);

f=figure;

f.Position(1)=f.Position(1)-400;
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*1.01; %0.785; %*0.4; %width
f.Position(4)=1170*0.47; %f %0.4 

t=tiledlayout(2,4);
%t=tiledlayout(2,3);
t.TileSpacing='tight';
%t.Padding='tight';
set(gcf,'color','w');

CT=flipud(cbrewer('div', 'RdBu', 30));

tile=(1:8); %6
tile=reshape(tile,4,2);

% labels=["a)" "b)" "c)" "d)"];
% labels_n=["e)" "f)" "g)"];

% loop through datasets 
for i=1:length(data_names)

    name=data_names(i);

    % file index
    index=file_i(3);

    % load file - noise added
    file=fullfile(d,filenames_d(index,:));
    data=load(file);
    % load file - no noise
    file_n=fullfile(d_n,filenames_n(index,:));
    data_n=load(file_n);

    % save time for 1st granule
    if isfield(data,'gn')

        day_1=char(data.gn.day(1));
        time_1=char(data.gn.airs_time(1)); 
        gn_1=char(data.gn.GN(1));
    end

    % % plot topography 
    % pcolor(ax(tile(i,1)),Topo_km.x, Topo_km.y, Topo_km.elev2'); shading flat
    % colormap(ax(tile(i,1)),"parula");
    % axis equal
    % xlim([-1500 2000]) 
    % ylim([-2000 2000]) 

    ax(tile(i,1))=nexttile(tile(i,1));
    pcolor(data.x, data.y, data.(name).Tp); shading flat 
    % change colormap 
    colormap(ax(tile(i,1)),CT);
    set(gca,'TickDir','out');
    clim([-3 3])
    % equal axes 
    axis equal
    xlim([-1500 2000]) % d11 %d19
    ylim([-2000 2000])             
    set(gca,'Color',[0.8 0.8 0.8])
    set(gca,'TickDir','out');
   
    title(names(i))
    label=labels(i,1);
    t_l=text(0,0,label);
    t_l.Units='normalized';
    t_l.Position=[-0.02 1.07];
    t_l.FontSize=10;

    set(gca,'Xticklabel',[]) 
    if i==1
       ylabel('Distance (km)') 
       % st=subtitle('AIRS Noise');
       % st.Rotation=90;
       % st.Position=[-2800 0];
       % st.FontWeight='bold';
       
       %l_n.Position=[-2750 1900];
    else
        set(gca,'Yticklabel',[])
        %l_n.Position=[-2200 1900];
    end
    ax(tile(i,1)).FontSize=9;


    % without added noise
    ax(tile(i,2))=nexttile(tile(i,2));
    if i==1 
        % plot topography
        pcolor(ax(tile(i,2)),Topo_km.x, Topo_km.y, Topo_km.elev2'); shading flat
        colormap(ax(tile(i,2)),map);
        set(gca,'Color','[0.1 0.55 0.8]')
        title('Topography')
        c=colorbar;
        c.FontSize=9;
        c.Label.String='Elevation (km)';
        c.Location='southoutside';
        c.Ticks = 0:6;
        c.TickLabels = num2cell(0:6) ;
        set(c,'TickDir','out');
        set(gca,'Color','[0.1 0.55 0.8]')
        hold on
        contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0]);
        hold off
    else
        pcolor(data.x, data.y, data_n.(name).Tp); shading flat 
        % change colormap 
        colormap(ax(tile(i,2)),CT);
        set(gca,'TickDir','out');
        clim([-3 3])
        title(strcat(names(i),' (no noise)'))
        set(gca,'Color',[0.8 0.8 0.8])
    end 
    % equal axes 
    axis equal
    xlim([-1500 2000]) % d11 %d19
    ylim([-2000 2000])                 
    set(gca,'TickDir','out');        
    xlabel('Distance (km)') 
    label=labels(i,2);
    t_l=text(0,0,label);
    t_l.Units='normalized';
    t_l.Position=[-0.02 1.07];
    t_l.FontSize=10;
    
    if i==1 
            ylabel('Distance (km)')
            % st=subtitle('No Noise');
            % 
            % st.Rotation=90;
            % st.Position=[-2800 0];
            % st.FontWeight='bold';
    else 
         set(gca,'Yticklabel',[])   
    end
    ax(tile(i,2)).FontSize=9;
    hold off
    



end

c=colorbar;
c.Layout.Tile='east';
c.Label.String='Temperature Perturbation (K)';
c.FontSize=9;
set(c,'TickDir','out');


% AIRS granule time
day=char(data.gn.day(end));
time=char(data.gn.airs_time(end));
GN=char(data.gn.GN(end));

t1=title(t, [time_1(1:2),':',time_1(4:5),' ',' - ', ...
    time(1:2),':', time(4:5),' ','UTC ',day(2),'^{th}',' November']);
t1.FontSize=12;

% save plot 
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps' ...
    '\example maps\tiled layout\colormap_models_Tp_',day,'_',time,'_2.png'],'Resolution',300);







