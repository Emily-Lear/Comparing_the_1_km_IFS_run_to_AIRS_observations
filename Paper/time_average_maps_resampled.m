%% Time average maps - IFS, AIRS and ERA5 

% Functions written by others: 
%%% bin2mat: Andrew Stevens, 3/10/2009
%%% cbrewer: Charles Robert, 06.12.2011

% granule stripes 

% Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night_prc\day_night_prc_70.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

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

%%

data_names=["ifs","s3d","era5","airs"];

% Directories 
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region noise no logs day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region noise no logs day_night\list_nighttime.mat';


day=load(f_day);
night=load(f_night);


data=struct('Day',day,'Night',night);
dn=["Day" "Night"];
field_names=["A","MF_z","MF_m"]; % "HW","VW","MF","AN"

% bins 
xbins=-4500:50:4500;
ybins=-2000:50:2000;
[yq, xq]=ndgrid(ybins,xbins);
sz_xq=size(xq);

store.xq=xq;
store.yq=yq;

store.N=zeros(sz_xq(1),sz_xq(2));

% day / night data
for i=1:2

    DN=dn(i);

    %loop through datasets
    for data_i=1:length(data_names)

        name=data_names{data_i};


        % make array of ones 
        point_num=ones(size(data.(DN).(name).A));
        % put 0s in locations of NaNs
        point_num(isnan(data.(DN).(name).A))=0;
   
        % bin ones in location of points 
        B1=bin2mat(data.(DN).x, data.(DN).y, point_num, xq,yq, '@sum');
        B1(isnan(B1))=0;
        store.(DN).(name).N=store.N+B1;


        % loop through fields 
        for j=1:length(field_names)
    
            field=field_names(j);
            
            % average angles 
            if strcmp(field,'AN')

                disp(1)

                % in radians 
                data.(DN).(name).AN=deg2rad(data.(DN).(name).AN);
    
                % bin data onto region
                store.(DN).(name).AN=bin2mat(data.(DN).x, data.(DN).y, data.(DN).(name).AN, xq,yq,'@circ_mean_omitnan');
                store.(DN).(name).AN=rad2deg(store.(DN).(name).AN);
    
                % change areas with no points back to NaNs
                store.(DN).(name).(field)(store.(DN).(name).N==0)=NaN;

           
            else

                data.(DN).(name).(field)(isnan(data.(DN).(name).(field)))=0;
        
                % bin data onto region
                store.(DN).(name).(field)=bin2mat(data.(DN).x, data.(DN).y, data.(DN).(name).(field), xq,yq,'@sum');
                % find mean
                store.(DN).(name).(field)=store.(DN).(name).(field)./...
                store.(DN).(name).N;
    
                % change areas with no points back to NaNs
                store.(DN).(name).(field)(store.(DN).(name).N==0)=NaN;

            end
            
        end 

    end

end

% % save binned data structure
save('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\time average maps\time_mean_day_night_MF.mat','-struct','store')


%% 70th prc amplitude cutoff time average map data

data_names=["ifs","s3d","era5","airs"];

% with amplitude cutoff
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat';

day=load(f_day);
night=load(f_night);


data=struct('Day',day,'Night',night);
dn=["Day" "Night"];
field_names="HW"; % "HW","VW","MF","AN"

% bins 
xbins=-4500:50:4500;
ybins=-2000:50:2000;
[yq, xq]=ndgrid(ybins,xbins);
sz_xq=size(xq);

store.xq=xq;
store.yq=yq;

store.N=zeros(sz_xq(1),sz_xq(2));

% day / night data
for i=1:2

    DN=dn(i);

    %loop through datasets
    for data_i=1:length(data_names)

        name=data_names{data_i};


        % make array of ones 
        point_num=ones(size(data.(DN).(name).A));
        % put 0s in locations of NaNs
        point_num(isnan(data.(DN).(name).A))=0;
   
        % bin ones in location of points 
        B1=bin2mat(data.(DN).x, data.(DN).y, point_num, xq,yq, '@sum');
        B1(isnan(B1))=0;
        store.(DN).(name).N=store.N+B1;


        % loop through fields 
        for j=1:length(field_names)
    
            field=field_names(j);
            
            % average angles 
            if strcmp(field,'AN')

                disp(1)

                % in radians 
                data.(DN).(name).AN=deg2rad(data.(DN).(name).AN);
    
                % bin data onto region
                store.(DN).(name).AN=bin2mat(data.(DN).x, data.(DN).y, data.(DN).(name).AN, xq,yq,'@circ_mean_omitnan');
                store.(DN).(name).AN=rad2deg(store.(DN).(name).AN);
    
                % change areas with no points back to NaNs
                store.(DN).(name).(field)(store.(DN).(name).N==0)=NaN;

           
            else

                data.(DN).(name).(field)(isnan(data.(DN).(name).(field)))=0;
        
                % bin data onto region
                store.(DN).(name).(field)=bin2mat(data.(DN).x, data.(DN).y, data.(DN).(name).(field), xq,yq,'@sum');
                % find mean
                store.(DN).(name).(field)=store.(DN).(name).(field)./...
                store.(DN).(name).N;
    
                % change areas with no points back to NaNs
                store.(DN).(name).(field)(store.(DN).(name).N==0)=NaN;

            end
            
        end 

    end

end

% % save binned data structure
save('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\time average maps\time_mean_day_night_70th_prc.mat','-struct','store')



%% Plot time average maps - noise

% filepath
f_data='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\time average maps\time_mean_day_night_MF.mat';
data_1=load(f_data);

% 70th prc
f_data_2='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\time average maps\time_mean_day_night_70th_prc.mat';
data_2=load(f_data_2);

dn=["Day" "Night"];
data_names=["airs","era5","ifs","s3d"];
data_names_2=["AIRS","ERA5 1","IFS 1","IFS 2"];
%fields=["A","HW","VW","MF","AN"];
fields=["A","HW","MF_z","MF_m"];
%field_names=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" ...
   % "Momentum Flux" "Angle"];
field_names=["Amplitude" "Horizontal Wavelength", "Zonal Momentum Flux" "Meridional Momentum Flux"];
%cb_labels=["Amplitude (K)" "Horizontal Wavelength (km)" "Vertical Wavelength (km)" ...
     %"log10( Momentum Flux )" ['Angle (',char(176),')']];
cb_labels=["Amplitude (K)" "Horizontal Wavelength (km)","Zonal Momentum Flux (mPa)" "Meridional Momentum Flux (mPa)"];
day_night=["Daytime" "Nighttime"];

labels=["(a)" "(b)" "(c)" "(d)"];
labels_2=["(e)" "(f)" "(g)" "(h)"];
labels_3=["(i)" "(j)" "(k)" "(l)"];

% cb_l_limits=[0 60 12 -3.5 -180];
% cb_u_limits=[2 430 35 -1.6 180];
% % Nighttime 
% cb_u_limits_n=[2.5 430 26 -1.6 180];

cb_l_limits=[0 60 -3.5 -3.5]; % MF: -2
cb_u_limits=[2 600 3.5 3.5]; % MF: 2
% Nighttime 
cb_u_limits_n=[2 430 3.5 3.5];

% [0 60 12 -5 -180]
% [3 490 45 -1.5 180]


% LonBox=[10 170];
% LatBox=[0 90];
% xRange=[-4700 4700];
% yRange=[-2200 2200];
% lon1=94;
% lat1=52;
% Loc1=[lat1 lon1];
% point_spacing=[30 30 1 1];
% % add topography 
% [Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);


% Colormaps
CT_1=cbrewer('seq', 'OrRd', 30);
% CT_2=cbrewer('seq', 'Blues', 30);
% CT_3=cbrewer('seq', 'Greens', 30);
% CT_4=cbrewer('seq', 'Purples', 30); 
% CT_5=cbrewer('div', 'BrBG', 30);
CT_2=cbrewer('seq', 'Blues', 30);
CT_3_a=flipud(cbrewer('div', 'BrBG',24)); % MF_z
CT_4_a=cbrewer('div', 'PRGn',24); % MF_m 
CT_3=ones(30,3);
CT_3(1:12,:)=CT_3_a(1:12,:);
CT_3(19:30,:)=CT_3_a(13:24,:);
CT_4=ones(30,3);
CT_4(1:12,:)=CT_4_a(1:12,:);
CT_4(19:30,:)=CT_4_a(13:24,:);

%CT_3(11:20,:)=1;

DN_names=["Daytime" "Nighttime"];

% CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5);
% CT_names = fieldnames(CT_s);
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4);
CT_names = fieldnames(CT_s);



% day/night 

for i=2

    DN=dn(i);
    DN_name=DN_names(i);

    % loop though fields
    for j=1:length(fields)

        cb_ll=cb_l_limits(j);

        if strcmp(DN,'Night')
            cb_ul=cb_u_limits_n(j);
        else
        cb_ul=cb_u_limits(j);
        end

        CT_name=CT_names{j};
        CT=CT_s.(CT_name);

        field=fields(j);
        if strcmp(field,'HW')
            data=data_2;
        else 
            data=data_1;
        end

%         disp(field)

        f=figure;
        %set(gcf,'WindowState','maximized')
        

        t=tiledlayout(2,2);
        t.TileSpacing='tight';
        t.Padding='compact';
        set(gcf,'color','w');
        f.Position(2)=f.Position(2)-100;
        f.Position(3)=f.Position(3)*2.1;
        f.Position(4)=f.Position(4)*1.4;

        % loop through datasets 
        for data_i=1:length(data_names)

            name=data_names(data_i);
            name_2=data_names_2(data_i);  

            %max_HW=max(data.(DN).(name).HW(:));

            % if strcmp(field,"MF")
            % data.(DN).(name).MF=log10(data.(DN).(name).MF);
            % end 


            if strcmp(field,"MF_z") || strcmp(field,"MF_m")
                
                data.(DN).(name).(field)=-(data.(DN).(name).(field)).*10^3;

                % log MF_z and MF_m
                MF=data.(DN).(name).(field);
                MF_sign=sign(MF);

                % log absolute values + 0.005
                log_MF=log10(abs(MF)+0.005);
                min_log_MF=min(log_MF(:));
                % make all values positive
                log_MF=log_MF+abs(log10(0.005)); 

                %log_MF(zero_i)=0;
                % change signs back 
                data.(DN).(name).(field)=(log_MF.*MF_sign);

                if strcmp(field,'MF_z')
                    min_MF_z=min(data.(DN).(name).MF_z(:));
                    max_MF_z=max(data.(DN).(name).MF_z(:));
                elseif strcmp(field,'MF_m')
                    min_MF_m=min(data.(DN).(name).MF_m(:));
                    max_MF_m=max(data.(DN).(name).MF_m(:));
                end

                % colorbar values
                cb_labels_2=[-1e2 -10 -1 -0.1 0 0.1 1 10 1e2];
                cb_ticks=[];

                for i=1:9

                    cbl=cb_labels_2(i);
                    cb_ticks(i)=(log10(abs(cbl)+0.005)+abs(log10(0.005)))*sign(cbl);

                end
                

            end
            if strcmp(name,"airs") && strcmp(field,'A')
                data.(DN).airs.(field)=(data.(DN).airs.(field))./2;
            end

            
%             disp(name)
%             disp(max(data.(DN).(name).(field)(:),[],'omitnan'))
%             disp(min(data.(DN).(name).(field)(:),[],'omitnan'))
            
            ax(data_i)=nexttile(data_i);
            set(gcf,'color','w');
            colormap(CT);
            pcolor(data.xq, data.yq,data.(DN).(name).(field)); 
            shading flat
            axis equal
            clim([cb_ll cb_ul])
            xlim([-4500 4500])
            ylim([-2000 2000])
            title(strjoin(name_2));
            set(gca,'Color',[0.8 0.8 0.8])
            set(gca,'TickDir','out');            
            if data_i==1 || data_i==3
                ylabel('Distance (km)')
            else
                 set(gca,'Yticklabel',[]) 
            end
            if data_i==3 || data_i==4
                xlabel('Distance (km)')
            else 
                 set(gca,'Xticklabel',[]) 
            end
            if strcmp(field,"MF_z")
                label=labels_2(data_i);
            elseif strcmp(field,"MF_m")
                 label=labels_3(data_i);
            else 
                label=labels(data_i);
            end
            t_l=text(0,0,label);
            t_l.Units='normalized';
            t_l.Position=[0.01 1.1];
            set(gca,'FontSize',16)  
            
        end

        c=colorbar;
        c.Label.String=cb_labels(j);
        set(c,'TickDir','out');
        set(c,'FontSize',16)
        c.Layout.Tile='east';
        %c.FontSize=12;
        if strcmp(field,"MF_z") || strcmp(field,"MF_m")
            c.Ticks=cb_ticks;
            c.TickLabels=cb_labels_2;
        end

       
        

        %title(t,strjoin([DN_name,' ',field_names(j)],''));
        title(t,[field_names(j)]);
        fontsize(16,"points")
        
        exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\time average maps\', ...
                    DN,' noise\time_mean_',field,'_paper.png'],''),'Resolution',300);
% 70th prc new
    end
end



