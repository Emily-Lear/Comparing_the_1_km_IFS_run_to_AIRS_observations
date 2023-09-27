% resampled data scatterplots - added AIRS noise


% 7x7 boxcar smoother

% momentum flux logged, less noise
% 
% % Percentiles
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night prc\day_night_prc_70.mat';
p=load(f_prc);

day_prc=2:2:20;
night_prc=3:2:21;

% data directory
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night';
% granule stripes day_night
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

data_names=["ifs","s3d","era5","airs"];
%data_names=fieldnames(data);

field_names=["A","HW","VW","AN","MF","T","MF_z","MF_m"];


%field_names=fieldnames(fields);

% % make stucture of indexes for each half day - 1-10 Nov, granule stripes
file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
    'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
    'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
    'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
    116:121, 'd21',123:127);


% % % make stucture of indexes for each half day (days 1-14)
% % missing data for d25
% file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
%     'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
%     'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
%     'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
%     116:121, 'd21',123:129, 'd22',130:136,'d23',137:142,'d24',143:147, ...
%     'd25',148,'d26',149:155,'d27',156:161,'d28',162:167,'d29',168:172 ...
%     ); 

half_days=fieldnames(file_indices);

% loop through half days of data
for i=2:length(half_days)

    %get file indexes 
    field_hd=half_days{i}; 
    file_i=file_indices.(field_hd);


    % make empty lists  for data 
    % loop through datasets 
    for data_j=1:length(data_names)
        name=data_names{data_j};
        % loop through field names 
        for field_j=1:length(field_names)
            field=field_names{field_j};

            list.(name).(field)=[];
        end 
    end 

    list.x=[];
    list.y=[];


    % loop through granules in half day 
    for j=1:length(file_i)
   
        % file index
        index=file_i(j);
        filename=filenames_d(index,:);
        % load file 
        f=fullfile(D,filename);
        data=load(f);



        % find indicies outside of x range
        x_i=find(data.x< -4500 | data.x>4500); %4500 %500 - 1500
        % find indicies outside of y range 
        y_i=find(data.y< -2000 | data.y>2000);

        % append x and y values to a list
        list.x=[list.x data.x(:)'];
        list.y=[list.y data.y(:)'];

        % AIRS granule time
        day=filename(9:10);
        time=filename(12:17);
    
        GN=filename(19:21);

        % save time for 1st granule
        if j==1

            filename_1=filenames_d(index,:);
            list.day_1=filename_1(9:10);
            list.time_1=filename_1(12:17);   

        end

        % loop through datasets 

        for data_i=1:length(data_names)

            name=data_names{data_i};

            %smooth granule with 7x7 filter 
            s.(name).A=smoothn(data.(name).A,[7 7]);


            if ismember(i,day_prc)
            amp_i=find(s.(name).A<p.(name).day);

            elseif ismember(i,night_prc)
            amp_i=find(s.(name).A<p.(name).night);
            end


            % loop through fields 
            for field_i=1:length(field_names)

                field=field_names{field_i};

                if strcmp(field,'MF')

                    % log momentum flux to base 10
                    data.(name).MF=log10((data.(name).MF).*10^3);

                end
        
                % remove points outside of region 
                data.(name).(field)(x_i)=NaN;
                data.(name).(field)(y_i)=NaN;

%               % remove points where amp <70th percentile K 
                data.(name).(field)(amp_i)=NaN;

                % concatenate data into list for half day
                list.(name).(field)=[list.(name).(field) data.(name).(field)(:)'];

                size(list.(name).(field));
         
            end
        end

    end
    
    % find indicies of nan values in all datasets for half days
    list.nan_i.i_nan=find(isnan(list.ifs.A));
    list.nan_i.s_nan=find(isnan(list.s3d.A));
    list.nan_i.e_nan=find(isnan(list.era5.A));
    list.nan_i.a_nan=find(isnan(list.airs.A));

    
    % % save lists for half day
    % % save structure 
    save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter\list_ST_',day,'_',time,'_',GN,'.mat'),'-struct', 'list');
%scatter region 70th filter
%scatter region noise 1_14

end

%% Plot scatter plots 
% save structure 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter';
%D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\box x_500-1500 km y_0-1000 km 1 K';
DI_d=dir(D);
filenames=char(DI_d.name);
filenames(1:2,:)=[];
sz_f=size(filenames);

field_names=["A","HW","VW","MF","AN"];

% Colormaps
CT_1=cbrewer('seq', 'OrRd', 30);
CT_2=cbrewer('seq', 'Blues', 30);
CT_3=cbrewer('seq', 'Greens', 30);
CT_4=cbrewer('seq', 'Purples', 30); 
CT_5=colormap(parula(30));
    
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5);
CT_names = fieldnames(CT_s);

% lower lim
xy_llimits=[0.5 60 12 -8 -180]; %0.61 %0.54 %0.78
xy_ulimits=[7 490 45 -1 180];

% xy_llimits=[-2 1.7 12 -8 -180]; %0.61 %0.54 %0.78
% xy_ulimits=[0 2.8 45 -1 180];

folders=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" "Momentum Flux" "Angle"];
corr_pval=["C" "P"];

% empty correlation arrays 
% day=2:2:20;
% night=1:2:21;

c_num=zeros(1,sz_f(1));

data1_data2=struct('C',c_num,'P',c_num);
fields=struct('A',data1_data2,'HW',data1_data2,'VW',data1_data2,'MF', ...
    data1_data2,'AN',data1_data2);

corr=struct('ifs1_ifs2',fields,'airs_era5',fields,'airs_ifs1', ...
    fields,'airs_ifs2',fields);


% loop through files in directory 
for i=1:length(filenames(:,1))

    % file index
    filename=filenames(i,:);
    % load file 
    f=fullfile(D,filename);
    list=load(f);

    time=filename(12:17);
    day_gn=filename(9:10);

    for j=1:length(field_names)

        field=field_names(j);
        CT_name=CT_names{j};
        CT=CT_s.(CT_name);
        xy_ll=xy_llimits(j);
        xy_ul=xy_ulimits(j);

%         if strcmp(field,"A") || strcmp(field,"HW") || strcmp(field,"MF")
% 
%             list.ifs.(field)=log10(list.ifs.(field));
%             list.era5.(field)=log10(list.era5.(field));
%             list.airs.(field)=log10(list.airs.(field));
%             list.s3d.(field)=log10(list.s3d.(field));
% 
%         end

        %figure
        
        % for the amplitude use the maximum value rounded up to the 
        % nearest integer
        if j==1
            % find max of wave property in all datasets for limits
            m_i=max(list.ifs.(field),[],'omitnan');
            m_e=max(list.era5.(field),[],'omitnan');
            m_a=max(list.airs.(field),[],'omitnan');
            m_s=max(list.s3d.(field),[],'omitnan');
            xy_ul=ceil(max([m_i m_e m_a m_s]));
        end

       % find min of wave property in all datasets for limits
        min_i=min(list.ifs.(field),[],'omitnan');
        min_e=min(list.era5.(field),[],'omitnan');
        min_a=min(list.airs.(field),[],'omitnan');
        min_s=min(list.s3d.(field),[],'omitnan');
        xy_llimit=floor(min([min_i min_e min_a min_s]));

        % find unique nan value indicies for all scatter plots 
        n_ifs_s3d=unique([list.nan_i.i_nan list.nan_i.s_nan]);
        n_airs_ifs=unique([list.nan_i.a_nan list.nan_i.i_nan]);
        n_airs_s3d=unique([list.nan_i.a_nan list.nan_i.s_nan]);
        n_airs_era5=unique([list.nan_i.a_nan list.nan_i.e_nan]);
        n_ifs_era5=unique([list.nan_i.i_nan list.nan_i.e_nan]);
        n_s3d_era5=unique([list.nan_i.s_nan list.nan_i.e_nan]);

        t=tiledlayout(2,3); %3
        t.TileSpacing='compact';
        t.Padding='tight';
        set(gcf,'color','w');
        
        % no. points in ifs 
        n=length(list.ifs.(field));

        if n>1

            % remove nan values
            ifs_s=list.ifs.(field)';
            ifs_s(n_ifs_s3d)=[];
            s3d_i=list.s3d.(field)';
            s3d_i(n_ifs_s3d)=[];

            n_is=length(s3d_i);

            ax(1)=nexttile(1);
            density_scatter(ifs_s,s3d_i)
            %dscatter(ifs_s,s3d_i,'smoothing',2500,'bins',[50,50],'marker','o');
            colormap(ax(1),CT);
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('IFS 1')
            ylabel('IFS 2')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            hold on 
            % line of best fit 
            p=polyfit(ifs_s, s3d_i,1);
            y=polyval(p,ifs_s);
            plot(ifs_s,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off
            
            % remove nan values
            airs_i=list.airs.(field)';
            airs_i(n_airs_ifs)=[];
            ifs_a=list.ifs.(field)';
            ifs_a(n_airs_ifs)=[];

            n_ai=length(airs_i);

            ax(4)=nexttile(4);
            density_scatter(airs_i,ifs_a)
            %dscatter(airs_i, ifs_a,'smoothing',2500,'bins',[50,50],'marker','o');
            colormap(ax(4),CT)
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('AIRS')
            ylabel('IFS 1')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            hold on
            % line of best fit 
            p=polyfit(airs_i,ifs_a,1);
            y=polyval(p,airs_i);
            plot(airs_i,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off

            % remove nan values
            airs_s=list.airs.(field)';
            airs_s(n_airs_s3d)=[];
            s3d_a=list.s3d.(field)';
            s3d_a(n_airs_s3d)=[];

            n_as=length(airs_s);
    
            ax(5)=nexttile(5);
            density_scatter(airs_s,s3d_a)
            %dscatter(airs_s, s3d_a,'smoothing',2500,'bins',[50,50],'marker','o');
            colormap(ax(5),CT)
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('AIRS')
            ylabel('IFS 2')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            % line of best fit
            hold on
            p=polyfit(airs_s, s3d_a,1);
            y=polyval(p,airs_s);
            plot(airs_s,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off

            % remove nan values
            airs_e=list.airs.(field)';
            airs_e(n_airs_era5)=[];
            era5_a=list.era5.(field)';
            era5_a(n_airs_era5)=[];

            n_ae=length(era5_a);
            
            ax(6)=nexttile(6);
            density_scatter(airs_e,era5_a)
            %dscatter(airs_e,era5_a,'smoothing',2500,'bins',[50,50],'marker','o');
            colormap(ax(6),CT)
%             c=colorbar;
%             c.Label.String='Density';
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('AIRS')
            ylabel('ERA5 1')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            % line of best fit 
            hold on
            p=polyfit(airs_e, era5_a,1);
            y=polyval(p,airs_e);
            plot(airs_e,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off

            % remove nan values
            ifs_e=list.ifs.(field)';
            ifs_e(n_ifs_era5)=[];
            era5_i=list.era5.(field)';
            era5_i(n_ifs_era5)=[];

            n_ie=length(ifs_e);
            
            ax(2)=nexttile(2);
            density_scatter(ifs_e,era5_i)
            %dscatter(ifs_e,era5_i,'smoothing',25,'bins',[200,200],'marker','o');
            colormap(ax(2),CT)
%             c=colorbar;
%             c.Label.String='Density';
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('IFS 1')
            ylabel('ERA5 1')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            % line of best fit 
            hold on
            p=polyfit(ifs_e, era5_i,1);
            y=polyval(p,ifs_e);
            plot(ifs_e,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off

            %remove nan values
            s3d_e=list.s3d.(field)';
            s3d_e(n_s3d_era5)=[];
            era5_s=list.era5.(field)';
            era5_s(n_s3d_era5)=[];

            n_se=length(s3d_e);
             
            ax(3)=nexttile(3);
            density_scatter(s3d_e,era5_s)
            %dscatter(s3d_e,era5_s,'smoothing',25,'bins',[200,200],'marker','o');
            colormap(ax(3),CT)
            c=colorbar;
            c.Label.String='Density';
            axis square
            set(gca,'Color',[0.92 0.92 0.92])
            xlabel('IFS 2')
            ylabel('ERA5 1')
            xlim([xy_ll xy_ul])
            ylim([xy_ll xy_ul])
            % line of best fit 
            hold on
            p=polyfit(s3d_e, era5_s,1);
            y=polyval(p,s3d_e);
            plot(s3d_e,y,'k');
            hold on
            % 1:1 line
            plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
            hold off
    
            t.Title.String=strjoin([field,' ',list.day_1,' ',list.time_1(1:2),':',list.time_1(4:5),' - ', ...
            ' ',day_gn,' ',time(1:2),':',time(4:5)],'');
    
            set(c,'TickDir','out');
            c.Layout.Tile='east';
    
            pause(1)
    
            folder=folders(j);
    
              % save plots 
              f=gcf;
%             exportgraphics(f,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\density scatter plots\box x_500-1500 km y_0-1000 km 1 K' ...
%             '/',folder,'/scatter_',day_gn,'_',time(1:2),'h',time(4:5),'m.jpg'],''));
              exportgraphics(f,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\density scatter plots\scatter region 70th filter' ...
                '/',folder,'/scatter_',day_gn,'_',time(1:2),'h',time(4:5),'m.png'],''),'Resolution',300);
       
            % find correlations between datasets 
            [C1,P1]=corrcoef(ifs_s,s3d_i);
            [C2,P2]=corrcoef(airs_i, ifs_a);
            [C3,P3]=corrcoef(airs_s, s3d_a);
            [C4,P4]=corrcoef(airs_e,era5_a);
            [C5,P5]=corrcoef(ifs_e,era5_i);
            [C6,P6]=corrcoef(s3d_e,era5_s);

            CP1=[C1(2),P1(2)];
            CP2=[C2(2),P2(2)];
            CP3=[C3(2),P3(2)];
            CP4=[C4(2),P4(2)];
            CP5=[C5(2),P5(2)];
            CP6=[C6(2),P6(2)];

            for k=1:2
                C_P=corr_pval(k);   
    
                corr.ifs1_ifs2.(field).(C_P)(i)=CP1(k);
                corr.airs_ifs1.(field).(C_P)(i)=CP2(k);
                corr.airs_ifs2.(field).(C_P)(i)=CP3(k);
                corr.airs_era5.(field).(C_P)(i)=CP4(k);  
                corr.ifs1_era5.(field).(C_P)(i)=CP5(k); 
                corr.ifs2_era5.(field).(C_P)(i)=CP6(k); 

            end

        end
    end
end 

% save correlation coeficients and P values 
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\corr 70th filter.mat'),'-struct', 'corr');

%save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\corr box x_500-1500 km y_0-1000 km 1 K.mat'),'-struct', 'corr');

 %'Resolution', 600

%% day and night data - 70th filter 
% loop through lists 
% save structure 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter';
%D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter';

DI_d=dir(D);
filenames=char(DI_d.name);
filenames(1:2,:)=[];
sz_f=size(filenames);

field_names=["A","HW","VW","MF","AN","T","MF_z" "MF_m"];

day_prc=1:2:20;
night_prc=2:2:20;

data_names=["ifs","s3d","era5","airs"];

% make empty lists  for data 
% loop through datasets 
for data_j=1:length(data_names)
    name=data_names{data_j};
    % loop through field names 
    for field_j=1:length(field_names)
        field=field_names{field_j};

        day.(name).(field)=[];
        night.(name).(field)=[];
    end 
end 

day.x=[];
day.y=[];
night.x=[];
night.y=[];

% loop through files in directory 
for i=1:length(filenames(:,1))

    % file index
    filename=filenames(i,:);
    % load file 
    f=fullfile(D,filename);
    list=load(f);

   if ismember(i,day_prc)
       day.x=[day.x list.x];
       day.y=[day.y list.y];
        
   else 
       night.x=[night.x list.x];
       night.y=[night.y list.y];
   end

    time=filename(12:17);
    day_gn=filename(9:10);

    for j=1:length(field_names)

        field=field_names(j);

        if ismember(i,day_prc)

            day.ifs.(field)=[day.ifs.(field) list.ifs.(field)];
            day.airs.(field)=[day.airs.(field) list.airs.(field)];
            day.era5.(field)=[day.era5.(field) list.era5.(field)];
            day.s3d.(field)=[day.s3d.(field) list.s3d.(field)];
        else

            night.ifs.(field)=[night.ifs.(field) list.ifs.(field)];
            night.airs.(field)=[night.airs.(field) list.airs.(field)];
            night.era5.(field)=[night.era5.(field) list.era5.(field)];
            night.s3d.(field)=[night.s3d.(field) list.s3d.(field)];

        end



    end

    

end

% find indicies of nan values in all datasets for daytime/nighttime
day.nan_i.i_nan=find(isnan(day.ifs.A));
day.nan_i.s_nan=find(isnan(day.s3d.A));
day.nan_i.e_nan=find(isnan(day.era5.A));
day.nan_i.a_nan=find(isnan(day.airs.A));

night.nan_i.i_nan=find(isnan(night.ifs.A));
night.nan_i.s_nan=find(isnan(night.s3d.A));
night.nan_i.e_nan=find(isnan(night.era5.A));
night.nan_i.a_nan=find(isnan(night.airs.A));

% save structure 
% save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat'),'-struct', 'day');
% save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat'),'-struct', 'night');
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat'),'-struct', 'day');
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat'),'-struct', 'night');


%% Plot scatter plots for daytime and nighttime data

% data files 
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat';

% f_stats='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\stats_day_night.mat';
% stats=load(f_stats);

day=load(f_day);
night=load(f_night);

data_names=["ifs" "s3d" "era5" "airs"];

% remove temperature
for k=1:length(data_names)
    name=data_names(k); 
    day.(name)=rmfield(day.(name),'T');
    night.(name)=rmfield(night.(name),'T');

end

clear day.ifs.T day.s3d.T day.era5.T day.airs.T ...
night.ifs.T night.s3d.T night.era5.T night.airs.T

dn=["Day" "Night"];
field_names=["A","HW","VW","MF","AN"];
corr_pval=["C" "P"];

data=struct('Day',day,'Night',night);

clear day night
% xy_llimits=[-2 1.7 12 -8 -180]; %0.61 %0.54 %0.78
% xy_ulimits=[0 2.8 45 -1 180];

% xy_llim_day=[0 60 12 -8 -180]; %0.61 %0.54 %0.78
% xy_ulim_day=[7 490 45 -1 180];
% 
% xy_lim_night=[0 60 12 -8 -180]; %0.61 %0.54 %0.78
% xy_ulim_night=[7 490 45 -1 180];

folders=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" "Momentum Flux" "Angle"];


% Colormaps
CT_1=cbrewer('seq', 'OrRd', 30);
CT_2=cbrewer('seq', 'Blues', 30);
CT_3=cbrewer('seq', 'Greens', 30);
CT_4=cbrewer('seq', 'Purples', 30); 
CT_5=colormap(parula(30));
    
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5);
CT_names = fieldnames(CT_s);
 % day / night
for i=1:2 

        DN=dn(i);

        if i==1
            xy_llimits=[0 60 10 -6 -180]; 
            xy_ulimits=[7 500 45 0 180];
        else 
            xy_llimits=[0 60 10 -6 -180]; 
            xy_ulimits=[15 500 46 0 180];    
        end

        % loop through fields 
        for j=1:length(field_names)

                field=field_names(j);
                CT_name=CT_names{j};
                CT=CT_s.(CT_name);
                xy_ll=xy_llimits(j);
                xy_ul=xy_ulimits(j);

                

                %figure 

                t=tiledlayout(2,3); %3
                t.TileSpacing='compact';
                t.Padding='tight';
                set(gcf,'color','w');
                 
                
                % no. points in ifs 
                n=length(data.(DN).ifs.(field));

                % find unique nan value indicies for all scatter plots 
                n_ifs_s3d=unique([data.(DN).nan_i.i_nan data.(DN).nan_i.s_nan]);
                n_airs_ifs=unique([data.(DN).nan_i.a_nan data.(DN).nan_i.i_nan]);
                n_airs_s3d=unique([data.(DN).nan_i.a_nan data.(DN).nan_i.s_nan]);
                n_airs_era5=unique([data.(DN).nan_i.a_nan data.(DN).nan_i.e_nan]);
                n_ifs_era5=unique([data.(DN).nan_i.i_nan data.(DN).nan_i.e_nan]);
                n_s3d_era5=unique([data.(DN).nan_i.s_nan data.(DN).nan_i.e_nan]);
           
                % remove nan values
                ifs_s=data.(DN).ifs.(field)';
                ifs_s(n_ifs_s3d)=[];
                s3d_i=data.(DN).s3d.(field)';
                s3d_i(n_ifs_s3d)=[];
        
                n_is=length(s3d_i);
 
%                 ax(1)=nexttile(1);
%                 density_scatter(ifs_s,s3d_i, 'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(ifs_s,s3d_i,'smoothing',2500,'bins',[50,50],'marker','o');
%                 colormap(ax(1),CT);
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('IFS 1')
%                 ylabel('IFS 2')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 hold on 
%                 % line of best fit 
%                 [p]=polyfit(ifs_s, s3d_i,1); % S, mu
%                 %[y,ifs1_ifs2.(DN).se.(field)]=polyval(p,ifs_s,S);
%                 y=polyval(p,ifs_s);
%                 plot(ifs_s,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 % text 
% %                 hold on
% %                 text(0,-2.5,strcat('n=',num2str(stats.ifs1_ifs2.(DN).num),' c=', ...
% %                     num2str(stats.ifs1_ifs2.(DN).corr.(field).C),' p=',num2str( ...
% %                     stats.ifs1_ifs2.(DN).corr.(field).P)),'Clipping','off')
%                 hold off

                % statistics 
                ifs1_ifs2.(DN).num=length(ifs_s);
                % points above and below the 1:1 line
                % below line x-y
                ifs1_ifs2.(DN).n_b.(field)=sum(ifs_s>s3d_i)/ifs1_ifs2.(DN).num;
                % above line y-x
                ifs1_ifs2.(DN).n_a.(field)=sum(s3d_i>ifs_s)/ifs1_ifs2.(DN).num;

                ifs1_ifs2.(DN).mean_i1.(field)=mean(ifs_s);
                ifs1_ifs2.(DN).std_i1.(field)=std(ifs_s);
                ifs1_ifs2.(DN).mean_i2.(field)=mean(s3d_i);
                ifs1_ifs2.(DN).std_i2.(field)=std(s3d_i);
                % y=mx+c
%                 ifs1_ifs2.(DN).m.(field)=p(1);
%                 ifs1_ifs2.(DN).c.(field)=p(2);

                % remove nan values
                ifs_e=data.(DN).ifs.(field)';
                ifs_e(n_ifs_era5)=[];
                era5_i=data.(DN).era5.(field)';
                era5_i(n_ifs_era5)=[];
        
                n_ie=length(ifs_e);
                
%                 ax(2)=nexttile(2);
%                 density_scatter(ifs_e,era5_i,'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(ifs_e,era5_i,'smoothing',25,'bins',[200,200],'marker','o');
%                 colormap(ax(2),CT)
%         %             c=colorbar;
%         %             c.Label.String='Density';
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('IFS 1')
%                 ylabel('ERA5 1')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p]=polyfit(ifs_e, era5_i,1);
%                 %[y,ifs1_era5.(DN).se.(field)]=polyval(p,ifs_e, S);
%                 y=polyval(p,ifs_e);
%                 plot(ifs_e,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 hold off

                % statistics 
                ifs1_era5.(DN).num=length(ifs_e);
                % points above and below the 1:1 line
                % below line x-y
                ifs1_era5.(DN).n_b.(field)=sum(ifs_e>era5_i)/ifs1_era5.(DN).num;
                % above line y-x
                ifs1_era5.(DN).n_a.(field)=sum(era5_i>ifs_e)/ifs1_era5.(DN).num;
                %disp(sum(ifs_e==era5_i));

                ifs1_era5.(DN).mean_i1.(field)=mean(ifs_e);
                ifs1_era5.(DN).std_i1.(field)=std(ifs_e);
                ifs1_era5.(DN).mean_e.(field)=mean(era5_i);
                ifs1_era5.(DN).std_e.(field)=std(era5_i);
                % y=mx+c
%                 ifs1_era5.(DN).m.(field)=p(1);
%                 ifs1_era5.(DN).c.(field)=p(2);


                %remove nan values
                s3d_e=data.(DN).s3d.(field)';
                s3d_e(n_s3d_era5)=[];
                era5_s=data.(DN).era5.(field)';
                era5_s(n_s3d_era5)=[];

                n_se=length(s3d_e);
                 
%                 ax(3)=nexttile(3);
%                 density_scatter(s3d_e,era5_s,'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(s3d_e,era5_s,'smoothing',25,'bins',[200,200],'marker','o');
%                 colormap(ax(3),CT)
% %                 c=colorbar;
% %                 c.Label.String='Density';
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('IFS 2')
%                 ylabel('ERA5 1')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p]=polyfit(s3d_e, era5_s,1);
%                 %[y,ifs2_era5.(DN).se.(field)]=polyval(p,s3d_e,S);
%                 y=polyval(p,s3d_e);
%                 plot(s3d_e,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 hold off

                % statistics 
                ifs2_era5.(DN).num=length(s3d_e);
                % points above and below the 1:1 line
                % below line x-y
                ifs2_era5.(DN).n_b.(field)=sum(s3d_e>era5_s)/ifs2_era5.(DN).num;
                % above line y-x
                ifs2_era5.(DN).n_a.(field)=sum(era5_s>s3d_e)/ifs2_era5.(DN).num;

                ifs2_era5.(DN).mean_i2.(field)=mean(s3d_e);
                ifs2_era5.(DN).std_i2.(field)=std(s3d_e);
                ifs2_era5.(DN).mean_e.(field)=mean(era5_s);
                ifs2_era5.(DN).std_e.(field)=std(era5_s);
                % y=mx+c
%                 ifs2_era5.(DN).m.(field)=p(1);
%                 ifs2_era5.(DN).c.(field)=p(2);

                
                % remove nan values
                airs_i=data.(DN).airs.(field)';
                airs_i(n_airs_ifs)=[];
                ifs_a=data.(DN).ifs.(field)';
                ifs_a(n_airs_ifs)=[];
        
                n_ai=length(airs_i);
        
%                 ax(4)=nexttile(4);
%                 density_scatter(airs_i,ifs_a, 'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(airs_i, ifs_a,'smoothing',2500,'bins',[50,50],'marker','o');
%                 colormap(ax(4),CT)
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('AIRS')
%                 ylabel('IFS 1')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 hold on
%                 % line of best fit 
%                 [p]=polyfit(airs_i,ifs_a,1); %,S
%                 %[y,a_i1_se]=polyval(p,airs_i,S);
%                 y=polyval(p,airs_i);
%                 plot(airs_i,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 hold off

                % statistics 
                airs_ifs1.(DN).num=length(airs_i);
                % points above and below the 1:1 line
                % below line x-y
                airs_ifs1.(DN).n_b.(field)=sum(airs_i>ifs_a)/airs_ifs1.(DN).num;
                % above line y-x
                airs_ifs1.(DN).n_a.(field)=sum(ifs_a>airs_i)/airs_ifs1.(DN).num;

                airs_ifs1.(DN).mean_a.(field)=mean(airs_i);
                airs_ifs1.(DN).std_a.(field)=std(airs_i);
                airs_ifs1.(DN).mean_i1.(field)=mean(ifs_a);
                airs_ifs1.(DN).std_i1.(field)=std(ifs_a);
                % y=mx+c
%                 airs_ifs1.(DN).m.(field)=p(1);
%                 airs_ifs1.(DN).c.(field)=p(2);

        
                % remove nan values
                airs_s=data.(DN).airs.(field)';
                airs_s(n_airs_s3d)=[];
                s3d_a=data.(DN).s3d.(field)';
                s3d_a(n_airs_s3d)=[];
        
                n_as=length(airs_s);
        
%                 ax(5)=nexttile(5);
%                 density_scatter(airs_s,s3d_a, 'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(airs_s, s3d_a,'smoothing',2500,'bins',[50,50],'marker','o');
%                 colormap(ax(5),CT)
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('AIRS')
%                 ylabel('IFS 2')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 % line of best fit
%                 hold on
%                 [p]=polyfit(airs_s, s3d_a,1);
%                 %[y,airs_ifs2.(DN).se.(field)]=polyval(p,airs_s,S);
%                 y=polyval(p,airs_s);
%                 plot(airs_s,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 hold off

                %statistics 
                airs_ifs2.(DN).num=length(airs_s);
                % points above and below the 1:1 line
                % below line x-y
                airs_ifs2.(DN).n_b.(field)=sum(airs_s>s3d_a)/airs_ifs2.(DN).num;
                % above line y-x
                airs_ifs2.(DN).n_a.(field)=sum(s3d_a>airs_s)/airs_ifs2.(DN).num;

                airs_ifs2.(DN).mean_a.(field)=mean(airs_s);
                airs_ifs2.(DN).std_a.(field)=std(airs_s);
                airs_ifs2.(DN).mean_i2.(field)=mean(s3d_a);
                airs_ifs2.(DN).std_i2.(field)=std(s3d_a);
%                 %y=mx+c
%                 airs_ifs2.(DN).m.(field)=p(1);
%                 airs_ifs2.(DN).c.(field)=p(2);

                
        
                % remove nan values
                airs_e=data.(DN).airs.(field)';
                airs_e(n_airs_era5)=[];
                era5_a=data.(DN).era5.(field)';
                era5_a(n_airs_era5)=[];
        
                n_ae=length(era5_a);
                
%                 ax(6)=nexttile(6);
%                 density_scatter(airs_e,era5_a, 'sz',2,'xy_range',[xy_ll xy_ul])
%                 %dscatter(airs_e,era5_a,'smoothing',2500,'bins',[50,50],'marker','o');
%                 colormap(ax(6),CT)
%                 c=colorbar;
%                 c.Label.String='Normalised density';
%                 axis square
%                 set(gca,'Color',[0.92 0.92 0.92])
%                 xlabel('AIRS')
%                 ylabel('ERA5 1')
%                 xlim([xy_ll xy_ul])
%                 ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p]=polyfit(airs_e, era5_a,1);
%                 %[y, airs_era5.(DN).se.(field)]=polyval(p,airs_e,S);
%                 y=polyval(p,airs_e);
%                 plot(airs_e,y,'k');
%                 hold on
%                 % 1:1 line
%                 plot(xy_ll:xy_ul,xy_ll:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
%                 hold off


                % statistics 
                airs_era5.(DN).num=length(airs_e);
                % points above and below the 1:1 line
                % below line x-y
                airs_era5.(DN).n_b.(field)=sum(airs_e>era5_a)/airs_era5.(DN).num;
                % above line y-x
                airs_era5.(DN).n_a.(field)=sum(era5_a>airs_e)/airs_era5.(DN).num;

                airs_era5.(DN).mean_a.(field)=mean(airs_e);
                airs_era5.(DN).std_a.(field)=std(airs_e);
                airs_era5.(DN).mean_e.(field)=mean(era5_a);
                airs_era5.(DN).std_e.(field)=std(era5_a);
                % y=mx+c
%                 airs_era5.(DN).m.(field)=p(1);
%                 airs_era5.(DN).c.(field)=p(2);


%                 t.Title.String=strjoin([field,'', DN]);
%         
%                 set(c,'TickDir','out');
%                 c.Layout.Tile='east';
% 
%                 folder=folders(j);
% 
%                 exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\density scatter plots\scatter region day_night\' ...
%                     'scatter_',field,'_', DN,'.png'],''),'Resolution',300);
%        
                % find correlations between datasets 
                [C1,P1]=corrcoef(ifs_s,s3d_i);
                [C2,P2]=corrcoef(airs_i, ifs_a);
                [C3,P3]=corrcoef(airs_s, s3d_a);
                [C4,P4]=corrcoef(airs_e,era5_a);
                [C5,P5]=corrcoef(ifs_e,era5_i);
                [C6,P6]=corrcoef(s3d_e,era5_s);
    
                CP1=[C1(2),P1(2)];
                CP2=[C2(2),P2(2)];
                CP3=[C3(2),P3(2)];
                CP4=[C4(2),P4(2)];
                CP5=[C5(2),P5(2)];
                CP6=[C6(2),P6(2)];
    
                for k=1:2
                    C_P=corr_pval(k);   
        
                    ifs1_ifs2.(DN).corr.(field).(C_P)=CP1(k);
                    airs_ifs1.(DN).corr.(field).(C_P)=CP2(k);
                    airs_ifs2.(DN).corr.(field).(C_P)=CP3(k);
                    airs_era5.(DN).corr.(field).(C_P)=CP4(k);  
                    ifs1_era5.(DN).corr.(field).(C_P)=CP5(k); 
                    ifs2_era5.(DN).corr.(field).(C_P)=CP6(k); 

                end
               
        end
end 

% % % make structure of statistics 
stats=struct("ifs1_ifs2",ifs1_ifs2,"airs_ifs1",airs_ifs1,"airs_ifs2", ...
airs_ifs2,"airs_era5",airs_era5,"ifs1_era5",ifs1_era5,"ifs2_era5", ifs2_era5);

% save statistics 
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\stats_day_night.mat'),'-struct', 'stats');
