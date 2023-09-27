% kernel distribution - resampled models and AIRS 

% lists of data - 70th percentile for each data set day/night
% 1st 10 days of November 2018
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat';
% unresampled
f_era5='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\not resampled\list_era5_airs_70th_prc_7.mat';
f_ifs='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\not resampled\list_ifs_airs_70th_prc_7.mat';

d_list=load(f_day);
n_list=load(f_night);
e_list=load(f_era5);
i_list=load(f_ifs);

fields=["A","HW","VW","MF"];
field_names=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" "Momentum Flux"];
x_labels=["Amplitude (K)", "Horizontal Wavelength (km)", "Vertical Wavelength (km)", "log_{10}(Momentum Flux (mPa))"];
bin_widths=[0.05 2 0.5 0.1];

% x_values_d=struct('A',0:0.125:4,'HW',60:20:500,'VW',10:2:45,'MF',-10.4:0.4:-0.4);
% x_values_n=struct('A',0:0.125:8,'HW',60:20:500,'VW',10:2:45,'MF',-10.4:0.4:-0.4);
x_values_d=struct('A',0:0.05:4,'HW',60:2:500,'VW',10:0.5:45,'MF',-10.4:0.1:-0.4);
x_values_n=struct('A',0:0.05:18,'HW',0:2:500,'VW',0:0.1:45,'MF',-3:0.1:3);
% A 0:0.05:8,VW 10:0.5:45, HW 500, MF -10.4:0.1:-0.4

% HW: 60:1:500

%% Day

f=figure;

t=tiledlayout(2,2);
t.TileSpacing='compact';
t.Padding='tight';
set(gcf,'color','w')

% loop through fields
for fieldi=1:length(fields)

    field=fields(fieldi);
    
%     % total point numbers
%     tc_ifs=sum(~isnan(d_list.ifs.(field)));
%     tc_s3d=sum(~isnan(d_list.s3d.(field)));
%     tc_era5=sum(~isnan(d_list.era5.(field)));
%     tc_airs=sum(~isnan(d_list.airs.(field)));

    pd_ifs=fitdist(d_list.ifs.(field)','Kernel');
    pd_s3d=fitdist(d_list.s3d.(field)','Kernel');
    pd_era5=fitdist(d_list.era5.(field)','Kernel');
    pd_airs=fitdist(d_list.airs.(field)','Kernel');

%     if fieldi==2
%         pd_ifs=fitdist(d_list.ifs.(field)','Kernel','Width',20);
%         pd_s3d=fitdist(d_list.s3d.(field)','Kernel','Width',20);
%         pd_era5=fitdist(d_list.era5.(field)','Kernel','Width',20);
%         pd_airs=fitdist(d_list.airs.(field)','Kernel','Width',20);
%     end

    x=x_values_d.(field);

    y_ifs=pdf(pd_ifs,x);
    y_s3d=pdf(pd_s3d,x);
    y_era5=pdf(pd_era5,x);
    y_airs=pdf(pd_airs,x);

    % maixmum values 
    max_ifs=max(y_ifs,[],'omitnan');
    max_s3d=max(y_s3d,[],'omitnan');
    max_era5=max(y_era5,[],'omitnan');
    max_airs=max(y_airs,[],'omitnan');

    ax(fieldi)=nexttile(fieldi);
    plot(x,y_ifs,'LineWidth',0.8); %./max_ifs
    title(field_names(fieldi));
    hold on 
    plot(x,y_s3d,'LineWidth',0.8); %./max_s3d
    hold on 
    plot(x,y_era5,'LineWidth',0.8); %./max_era5
    hold on 
    plot(x,y_airs,'LineWidth',0.8); %./max_airs
    ylim([0 0.5])
    ax(fieldi).FontSize=10;
    hold off

%     box off
%     set(gca,'TickDir','out');
    xlabel(string(x_labels(fieldi)))
    
end

f.Position(1)=f.Position(1)-50;
f.Position(2)=f.Position(2)-80;
f.Position(3:4)=f.Position(3:4)*1.2;

title(t,'Daytime')
leg=legend('IFS 1','IFS 2','ERA5 1','AIRS','Orientation','horizontal');
leg.Layout.Tile='South';


% save figure 
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\kernel distribution\region day_night 70th prc\' ...
    'kbf_day.png'],'Resolution',300);


%% Night 

labels=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)"];
labels=reshape(labels,2,4)';

f=figure;

t=tiledlayout(4,2);
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w')
t.Position=[0.12 0.13 0.85 0.83]; %0.1

% x axis limits
x_ll=[0 50 5 -3];
x_ul=[18 500 45 3];
y_ll=[-8 0 0 0];
y_ul=[0 4 8 13];

TL=reshape(1:8,2,4)';

% resampled
% loop through fields 
for fieldi=1:length(fields)


    field=fields(fieldi);
    bw=bin_widths(fieldi);

%     % total point numbers
%     tc_ifs=sum(~isnan(n_list.ifs.field));
%     tc_s3d=sum(~isnan(n_list.s3d.field));
%     tc_era5=sum(~isnan(n_list.era5.field));
%     tc_airs=sum(~isnan(n_list.airs.field));

    % % log amplitude before finding kdf  
    % n_list.ifs.A=log10(n_list.ifs.A);
    % n_list.s3d.A=log10(n_list.s3d.A);
    % n_list.era5.A=log10(n_list.era5.A);
    % n_list.airs.A=log10(n_list.airs.A);

   
    % remove nan values
    n_list.ifs.(field)(isnan(n_list.ifs.(field)))=[];
    n_list.s3d.(field)(isnan(n_list.s3d.(field)))=[];
    n_list.era5.(field)(isnan(n_list.era5.(field)))=[];
    n_list.airs.(field)(isnan(n_list.airs.(field)))=[];

    pd_ifs=fitdist(n_list.ifs.(field)','Kernel');
    pd_s3d=fitdist(n_list.s3d.(field)','Kernel');
    pd_era5=fitdist(n_list.era5.(field)','Kernel');
    pd_airs=fitdist(n_list.airs.(field)','Kernel');



    x=x_values_n.(field);

    y_ifs=pdf(pd_ifs,x)*bw;
    y_s3d=pdf(pd_s3d,x)*bw;
    y_era5=pdf(pd_era5,x)*bw;
    y_airs=pdf(pd_airs,x)*bw;

    if fieldi==1
        y_ifs=log10(y_ifs);
        y_s3d=log10(y_s3d);
        y_era5=log10(y_era5);
        y_airs=log10(y_airs);
    else
        y_ifs=y_ifs*(10^(2));
        y_s3d=y_s3d*(10^(2));
        y_era5=y_era5*(10^(2));
        y_airs=y_airs*(10^(2));
    end

    % maximum values 
    max_ifs1=max(y_ifs,[],'omitnan');
    max_s3d2=max(y_s3d,[],'omitnan');
    max_era51=max(y_era5,[],'omitnan');
    max_airs=max(y_airs,[],'omitnan'); 

    ax_i=TL(fieldi,2);
    
    ax(ax_i)=nexttile(ax_i);
    h(1)=plot(x,y_ifs,'LineWidth',0.8); %./max_ifs
    
    hold on 
    h(2)=plot(x,y_s3d,'LineWidth',0.8); %./max_s3d
    hold on 
    h(3)=plot(x,y_era5,'LineWidth',0.8); %./max_era5
    hold on 
    h(4)=plot(x,y_airs,'LineWidth',0.8); %./max_airs
    xlim([x_ll(fieldi) x_ul(fieldi)])
    ylim([y_ll(fieldi) y_ul(fieldi)])
    set(gca,'YTickLabel',[]);
    if fieldi==1
        title('Resampled as AIRS');
    end

    %     ylabel('log_{10}(Probability)')
    % else
    %     ylabel('Probability \times 10^{-2}')
    % end  
    ax(ax_i).FontSize=9;
    set(gca,'TickDir','out');
    label=labels(fieldi,2);
    t_l=text(0,0,label);
    t_l.Units='normalized';
    t_l.Position=[0.92 0.87];
    hold off
        
%     box off
%     set(gca,'TickDir','out')
    xlabel(string(x_labels(fieldi)))
    
    
end

field_names_2=["Amplitude (no resampling)" "Horizontal Wavelength (no resampling)" "Vertical Wavelength (no resampling)" "Momentum Flux (no resampling)"];

% Before resampling
% loop through fields 
for fieldi=1:length(fields)


    field=fields(fieldi);
    bw=bin_widths(fieldi);

%     % total point numbers
%     tc_ifs=sum(~isnan(n_list.ifs.field));
%     tc_s3d=sum(~isnan(n_list.s3d.field));
%     tc_era5=sum(~isnan(n_list.era5.field));
%     tc_airs=sum(~isnan(n_list.airs.field));

    % % log amplitude before finding kdf  
    % i_list.A.night=log10(i_list.A.night);
    % e_list.A.night=log10(e_list.A.night);
    % n_list.airs.A=log10(n_list.airs.A);
    
    % remove nan values
    n_list.airs.(field)(isnan(n_list.airs.(field)))=[];
    e_list.(field).night(isnan(e_list.(field).night))=[];
    i_list.(field).night(isnan(i_list.(field).night))=[];

    pd_ifs=fitdist(i_list.(field).night','Kernel');
    pd_era5=fitdist(e_list.(field).night','Kernel');
    pd_airs=fitdist(n_list.airs.(field)','Kernel');

    x=x_values_n.(field);

    y_ifs=pdf(pd_ifs,x)*bw;
    y_era5=pdf(pd_era5,x)*bw;
    y_airs=pdf(pd_airs,x)*bw;

     if fieldi==1
        y_ifs=log10(y_ifs);
        y_era5=log10(y_era5);
        y_airs=log10(y_airs);
     else
        y_ifs=y_ifs*(10^(2));
        y_era5=y_era5*(10^(2));
        y_airs=y_airs*(10^(2));
     end

    % maximum values 
    max_ifs=max(y_ifs,[],'omitnan');
    max_era5=max(y_era5,[],'omitnan');
    max_airs=max(y_airs,[],'omitnan'); 

    ax_i=TL(fieldi,1);
    
    ax(ax_i)=nexttile(ax_i);
    h(5)=plot(x,y_ifs,'LineWidth',0.8,'Color',"#4DBEEE"); %./max_ifs  
    hold on 
    h(6)=plot(x,y_era5,'LineWidth',0.8,'Color',"#77AC30"); %./max_era5
    hold on 
    plot(x,y_airs,'LineWidth',0.8,'Color',"#7E2F8E"); %./max_airs
    xlim([x_ll(fieldi) x_ul(fieldi)])
    ylim([y_ll(fieldi) y_ul(fieldi)])
    if fieldi==1
        ylabel('log_{10}(Probability)')
        title('Not resampled'); %field_names_2(fieldi)
    else
        ylabel('Probability \times 10^{-2}')
    end
    st=subtitle(ax(ax_i),field_names(fieldi),'FontWeight','bold');
    st.Rotation=90;  
    st.Units='normalized';
    st.Position=[-0.16 0.5];
    label=labels(fieldi,1);
    t_l=text(0,0,label);
    t_l.Units='normalized';
    t_l.Position=[0.92 0.87];
    ax(ax_i).FontSize=9;
    set(gca,'TickDir','out');
    hold off
%     set(gca,'TickDir','out')
    xlabel(string(x_labels(fieldi)))
    
    
end


f.Position(1)=f.Position(1)-50;
f.Position(2)=f.Position(2)-200;
f.Position(3)=f.Position(3)*1.23;
f.Position(4)=f.Position(4)*1.6;

%title(t,'Nighttime')
leg=legend(h,'IFS 1','IFS 2','ERA5 1','AIRS','IFS','ERA5','Orientation','horizontal');
leg.Layout.Tile='South';

% save figure 
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\kernel distribution\region day_night 70th prc\' ...
    'kdf_night_log_A.png'],'Resolution',300);

