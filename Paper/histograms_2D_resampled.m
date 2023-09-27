%% 2D histograms - day/night plots

% data files 
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat';

f_stats='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\stats_day_night.mat';
stats=load(f_stats);

day=load(f_day);
night=load(f_night);

dn=["Day" "Night"];
field_names=["A","HW","VW","MF","AN"];
field_names_2=["A","HW","VW","log 10(MF)","AN"];
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
            %logged
%               xy_llimits=[-1.85 1.75 10 -6 -180]; 
%               xy_ulimits=[0.85 2.75 46 0 180];
        else 
            xy_llimits=[0 60 10 -6 -180]; 
            xy_ulimits=[15 500 46 0 180]; 
%             xy_llimits=[-1.85 1.75 10 -6 -180]; 
%             xy_ulimits=[1.25 2.72 46 0 180];    
        end

        % loop through fields 
        for j=1:length(field_names)

                field=field_names(j);
                CT_name=CT_names{j};
                CT=CT_s.(CT_name);
                xy_ll=xy_llimits(j);
                xy_ul=xy_ulimits(j);


                

%                 ifs_max=max(data.(DN).ifs.(field))
%                 ifs_min=min(data.(DN).ifs.(field))
%                 s3d_max=max(data.(DN).s3d.(field))
%                 s3d_min=min(data.(DN).s3d.(field))
%                 airs_max=max(data.(DN).airs.(field))
%                 airs_min=min(data.(DN).airs.(field))
%                 era5_max=max(data.(DN).era5.(field))
%                 era5_min=min(data.(DN).era5.(field))


                figure

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
 
                ax(1)=nexttile(1);
                density_scatter(ifs_s,s3d_i, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(ifs_s,s3d_i,'smoothing',2500,'bins',[50,50],'marker','o');
                colormap(ax(1),CT);
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('IFS 1')
                ylabel('IFS 2')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
                hold on 
%                 % line of best fit 
%                 [p,S]=polyfit(ifs_s, s3d_i,1); % S, mu
%                 [y,ifs1_ifs2.(DN).se.(field)]=polyval(p,ifs_s,S);
%                 %y=polyval(p,ifs_s);
%                 plot(ifs_s,y,'k');
%                 hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                % text 
%                 hold on
%                 % x axes length
%                 x_length=xy_ul-xy_ll;
%                 text((-0.05*x_length),(-0.45*x_length),strcat('n=',num2str(stats.ifs1_ifs2.(DN).num),' n\_a=', ...
%                 num2str(stats.ifs1_ifs2.(DN).n_a), '\newline', 'n\_b=',num2str(stats.ifs1_ifs2.(DN).n_b), ...
%                 ' c=', num2str(stats.ifs1_ifs2.(DN).corr.(field).C),' p=',...
%                 num2str(stats.ifs1_ifs2.(DN).corr.(field).P)),'Clipping','off','FontSize',8);
                hold off

                % remove nan values
                ifs_e=data.(DN).ifs.(field)';
                ifs_e(n_ifs_era5)=[];
                era5_i=data.(DN).era5.(field)';
                era5_i(n_ifs_era5)=[];
        
                n_ie=length(ifs_e);
                
                ax(2)=nexttile(2);
                density_scatter(ifs_e,era5_i,'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(ifs_e,era5_i,'smoothing',25,'bins',[200,200],'marker','o');
                colormap(ax(2),CT)
        %             c=colorbar;
        %             c.Label.String='Density';
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('IFS 1')
                ylabel('ERA5 1')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p,S]=polyfit(ifs_e, era5_i,1);
%                 [y,ifs1_era5.(DN).se.(field)]=polyval(p,ifs_e, S);
%                 %y=polyval(p,ifs_e);
%                 plot(ifs_e,y,'k');
                hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                hold off


                %remove nan values
                s3d_e=data.(DN).s3d.(field)';
                s3d_e(n_s3d_era5)=[];
                era5_s=data.(DN).era5.(field)';
                era5_s(n_s3d_era5)=[];

                n_se=length(s3d_e);
                 
                ax(3)=nexttile(3);
                density_scatter(s3d_e,era5_s,'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(s3d_e,era5_s,'smoothing',25,'bins',[200,200],'marker','o');
                colormap(ax(3),CT)
%                 c=colorbar;
%                 c.Label.String='Density';
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('IFS 2')
                ylabel('ERA5 1')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p,S]=polyfit(s3d_e, era5_s,1);
%                 [y,ifs2_era5.(DN).se.(field)]=polyval(p,s3d_e,S);
%                 %y=polyval(p,s3d_e);
%                 plot(s3d_e,y,'k');
                hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                hold off


                
                % remove nan values
                airs_i=data.(DN).airs.(field)';
                airs_i(n_airs_ifs)=[];
                ifs_a=data.(DN).ifs.(field)';
                ifs_a(n_airs_ifs)=[];
        
                n_ai=length(airs_i);
        
                ax(4)=nexttile(4);
                density_scatter(airs_i,ifs_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');                
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(airs_i, ifs_a,'smoothing',2500,'bins',[50,50],'marker','o');
                colormap(ax(4),CT)
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('AIRS')
                ylabel('IFS 1')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
                hold on
%                 % line of best fit 
%                 [p,S]=polyfit(airs_i,ifs_a,1);
%                 [y,airs_ifs1.(DN).se.(field)]=polyval(p,airs_i,S);
%                 %y=polyval(p,airs_i);
%                 plot(airs_i,y,'k');
%                 hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                hold off

        
                % remove nan values
                airs_s=data.(DN).airs.(field)';
                airs_s(n_airs_s3d)=[];
                s3d_a=data.(DN).s3d.(field)';
                s3d_a(n_airs_s3d)=[];
        
                n_as=length(airs_s);
        
                ax(5)=nexttile(5);
                density_scatter(airs_s,s3d_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(airs_s, s3d_a,'smoothing',2500,'bins',[50,50],'marker','o');
                colormap(ax(5),CT)
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('AIRS')
                ylabel('IFS 2')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
%                 % line of best fit
%                 hold on
%                 [p,S]=polyfit(airs_s, s3d_a,1);
%                 [y,airs_ifs2.(DN).se.(field)]=polyval(p,airs_s,S);
%                 %y=polyval(p,airs_s);
%                 plot(airs_s,y,'k');
                hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                hold off
                
        
                % remove nan values
                airs_e=data.(DN).airs.(field)';
                airs_e(n_airs_era5)=[];
                era5_a=data.(DN).era5.(field)';
                era5_a(n_airs_era5)=[];
        
                n_ae=length(era5_a);
                
                ax(6)=nexttile(6);  
                density_scatter(airs_e,era5_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
%                 % log normalised bin counts 
%                 total_counts=sum(h.BinCounts(:));
%                 counts=h.BinCounts./total_counts;
%                 counts(counts==0)=NaN;               
%                 counts_log=log10(counts);                  
%                 pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
                shading flat; box off; grid off
                %dscatter(airs_e,era5_a,'smoothing',2500,'bins',[50,50],'marker','o');
                colormap(ax(6),CT)
                c=colorbar;
                c.FontSize=8;
                if strcmp(field,'A') 
                c.Label.String='log10( Normalised density )';
                else
                    c.Label.String='Noramilsed density';
                end
                axis square
                box off
                set(gca,'TickDir','out');
                set(gca,'Color',[0.92 0.92 0.92])
                %set(gca,'ColorScale','log')
                xlabel('AIRS')
                ylabel('ERA5 1')
                xlim([xy_ll xy_ul])
                ylim([xy_ll xy_ul])
%                 % line of best fit 
%                 hold on
%                 [p,S]=polyfit(airs_e, era5_a,1);
%                 [y, airs_era5.(DN).se.(field)]=polyval(p,airs_e,S);
%                 %y=polyval(p,airs_e);
%                 plot(airs_e,y,'k');
                hold on
                % 1:1 line
                plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
                hold off

%                 % log A and HW 
%                 if strcmp(field,'A') 
% 
%                     set(ax(1),'ColorScale','log')
%                     set(ax(2),'ColorScale','log')
%                     set(ax(3),'ColorScale','log')
%                     set(ax(4),'ColorScale','log')
%                     set(ax(5),'ColorScale','log')
%                     set(ax(6),'ColorScale','log')
% %                     
%                 disp(1)
%                 end

                field_name=field_names_2(j);
                t.Title.String=strjoin([field_name,'', DN]);
        
                set(c,'TickDir','out');
                c.Layout.Tile='east';

                folder=folders(j);

                % exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\2D histograms\region day_night 70th prc 50x50 bins\' ...
                %     '2D_hist_',field,'_', DN,'.png'],''),'Resolution',300);
               
        end

        %exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\2D histograms\day_night plots 70th prc\' ...
                %     '2D_hist_', DN,'.png'],''),'Resolution',300);
end 

%% Plot day/night figures

% data files 
f_day='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_daytime.mat';
f_night='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region 70th filter day_night\list_nighttime.mat';

f_stats='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\correlations\stats_day_night.mat';
stats=load(f_stats);

day=load(f_day);
night=load(f_night);

dn=["Day" "Night"];
day_night=["Daytime" "Nighttime"];
field_names=["A","HW","VW","MF"]; %"AN"
field_names_2=["A (K)","HW (km)","VW (km)","log_{10}(MF (mPa))"];% "AN"
corr_pval=["C" "P"];

labels=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)" "(m)" "(n)" "(o)" "(p)" ...
    "(q)" "(r)" "(s)" "(t)" "(u)" "(v)" "(w)" "(x)"];
labels=reshape(labels,4,6);

data=struct('Day',day,'Night',night);

clear day night


folders=["Amplitude" "Horizontal Wavelength" "Vertical Wavelength" "Momentum Flux"];
% "Angle"

% tiled layout indicies 
TL=reshape(1:24,4,6);


% Colormaps
CT_1=cbrewer('seq', 'OrRd', 30);
CT_2=cbrewer('seq', 'Blues', 30);
CT_3=cbrewer('seq', 'Greens', 30);
CT_4=cbrewer('seq', 'Purples', 30); 
CT_5=colormap(parula(30));
    
CT_s = struct('CT_1',CT_1,'CT_2',CT_2,'CT_3',CT_3,'CT_4',CT_4,'CT_5',CT_5);
CT_names = fieldnames(CT_s);
 % day / night
for i=2 

    DN=dn(i);

    if i==1
        xy_llimits=[0 60 10 -5.2 -180]; 
        xy_ulimits=[7 500 45 0 180];
        c_llimits=[0 0 0 0];
        
        % logged
              % xy_llimits=[-1.85 60 10 -6 -180]; 
              % xy_ulimits=[0.85 500 45 0 180];
    else 
        xy_llimits=[0 60 10 -3 -180]; %-5.2
        xy_ulimits=[15 500 46 2.6 180]; %-0.2
        c_llimits=[0 0 0 0];
        c_ulimits=[0.018 0.010 0.014 0.014];
            % xy_llimits=[-1.85 60 10 -6 -180]; 
            % xy_ulimits=[1.25 500 46 0 180];    
    end

    f=figure;
    

    t=tiledlayout(6,4); % 6 rows, 4 columns
    t.TileSpacing='tight';
    t.Padding='compact';
    set(gcf,'color','w');
    %set(gcf,'Units','inches');
    f.Position(2)=f.Position(2)-400;
    f.Position(3)=830*0.54; %*0.4; %width
    f.Position(4)=1170*0.59; %f.Position(4)*11.7*0.4; %height



    % loop through fields 
    for j=1:length(field_names)
        
        % % log amplitudes 
        % if j==1
        %     data.(DN).ifs.A=log10(data.(DN).ifs.A);
        %     data.(DN).s3d.A=log10(data.(DN).s3d.A);
        %     data.(DN).airs.A=log10(data.(DN).airs.A);
        %     data.(DN).era5.A=log10(data.(DN).era5.A);
        % end

        field=field_names(j);
        CT_name=CT_names{j};
        CT=CT_s.(CT_name);
        xy_ll=xy_llimits(j);
        xy_ul=xy_ulimits(j);
        c_ll=c_llimits(j);
        c_ul=c_ulimits(j);

        

        % tiled layout column
        col=squeeze(TL(j,:));
        % labels column
        l_col=squeeze(labels(j,:));
         
        
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

        ax(col(1))=nexttile(col(1));
        if strcmp(field,'A') 
            ds=density_scatter(ifs_s,s3d_i, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            density_scatter(ifs_s,s3d_i, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
        end
        % min_1=min(ds.Values(:))
        % max_1=max(ds.Values(:))
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        shading flat; box off; grid off
        ax(col(1)).FontSize=6.5;
        %dscatter(ifs_s,s3d_i,'smoothing',2500,'bins',[50,50],'marker','o');
        colormap(ax(col(1)),CT);        
        ax(col(1)).CLim=[c_ll c_ul];
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88])
        %set(gca,'ColorScale','log')
        xlabel('IFS 1')
        if j==1
            ylabel('IFS 2')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        set(gca,'Xticklabel',[])
        title(field_names_2(j)) %field_names_2
        label=l_col(1);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on 
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off
        

        % remove nan values
        ifs_e=data.(DN).ifs.(field)';
        ifs_e(n_ifs_era5)=[];
        era5_i=data.(DN).era5.(field)';
        era5_i(n_ifs_era5)=[];

        n_ie=length(ifs_e);
        
        ax(col(2))=nexttile(col(2));
        if strcmp(field,'A') 
            ds=density_scatter(ifs_e,era5_i, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            density_scatter(ifs_e,era5_i,'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
        end
        % min_2=min(ds.Values(:))
        % max_2=max(ds.Values(:))
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        shading flat; box off; grid off
        ax(col(2)).FontSize=6.5;
        %dscatter(ifs_e,era5_i,'smoothing',25,'bins',[200,200],'marker','o');
        colormap(ax(col(2)),CT)
        clim([c_ll c_ul])       
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88])
        %set(gca,'ColorScale','log')
        xlabel('IFS 1')
        if j==1
            ylabel('ERA5 1')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        set(gca,'Xticklabel',[])
        label=l_col(2);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off
        

        %remove nan values
        s3d_e=data.(DN).s3d.(field)';
        s3d_e(n_s3d_era5)=[];
        era5_s=data.(DN).era5.(field)';
        era5_s(n_s3d_era5)=[];

        n_se=length(s3d_e);
         
        ax(col(3))=nexttile(col(3));
        if strcmp(field,'A') 
            ds=density_scatter(s3d_e,era5_s, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            density_scatter(s3d_e,era5_s,'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
        end
        % min_3=min(ds.Values(:))
        % max_3=max(ds.Values(:)) 
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        ax(col(3)).FontSize=6.5;
        shading flat; box off; grid off
        %dscatter(s3d_e,era5_s,'smoothing',25,'bins',[200,200],'marker','o');
        colormap(ax(col(3)),CT)
        clim([c_ll c_ul])
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88])
        %set(gca,'ColorScale','log')
        xlabel('IFS 2')
        if j==1
            ylabel('ERA5 1')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        set(gca,'Xticklabel',[])
        label=l_col(3);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off


        
        % remove nan values
        airs_i=data.(DN).airs.(field)';
        airs_i(n_airs_ifs)=[];
        ifs_a=data.(DN).ifs.(field)';
        ifs_a(n_airs_ifs)=[];

        n_ai=length(airs_i);

        ax(col(4))=nexttile(col(4));
        if strcmp(field,'A') 
            ds=density_scatter(airs_i,ifs_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            density_scatter(airs_i,ifs_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');                
        end
        % min_4=min(ds.Values(:))
        % max_4=max(ds.Values(:))
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        shading flat; box off; grid off
        ax(col(4)).FontSize=6.5;
        %dscatter(airs_i, ifs_a,'smoothing',2500,'bins',[50,50],'marker','o');
        colormap(ax(col(4)),CT)
        clim([c_ll c_ul])
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88])
        %set(gca,'ColorScale','log')
        xlabel('AIRS')
        if j==1
            ylabel('IFS 1')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        set(gca,'Xticklabel',[])
        label=l_col(4);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off
        


        % remove nan values
        airs_s=data.(DN).airs.(field)';
        airs_s(n_airs_s3d)=[];
        s3d_a=data.(DN).s3d.(field)';
        s3d_a(n_airs_s3d)=[];

        n_as=length(airs_s);

        ax(col(5))=nexttile(col(5));
        if strcmp(field,'A') 
            ds=density_scatter(airs_s,s3d_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            density_scatter(airs_s,s3d_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
        end
        % min_5=min(ds.Values(:))
        % max_5=max(ds.Values(:))
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        shading flat; box off; grid off
        ax(col(5)).FontSize=6.5;
        %dscatter(airs_s, s3d_a,'smoothing',2500,'bins',[50,50],'marker','o');
        colormap(ax(col(5)),CT)
        clim([c_ll c_ul])
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88]) %0.92 0.92 0.92
        %set(gca,'ColorScale','log')
        xlabel('AIRS')
        if j==1
            ylabel('IFS 2')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        set(gca,'Xticklabel',[])
        label=l_col(5);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off
        
        

        % remove nan values
        airs_e=data.(DN).airs.(field)';
        airs_e(n_airs_era5)=[];
        era5_a=data.(DN).era5.(field)';
        era5_a(n_airs_era5)=[];

        n_ae=length(era5_a);
        
        ax(col(6))=nexttile(col(6)); 
        if strcmp(field,'A') 
            ds=density_scatter(airs_e,era5_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'log_10_hist');
            max(ds.Values(:));
        else 
            ds=density_scatter(airs_e,era5_a, 'sz',2,'xy_range',[xy_ll xy_ul],'num_bins',[50 50],'2D_hist');
        end
        % min_6=min(ds.Values(:))
        % max_6=max(ds.Values(:))
        % % log normalised bin counts 
        % total_counts=sum(h.BinCounts(:));
        % counts=h.BinCounts./total_counts;
        % counts(counts==0)=NaN;               
        % counts_log=log10(counts);                  
        % pcolor(h.XBinEdges(1:end-1),h.YBinEdges(1:end-1),counts_log')
        shading flat; box off; grid off
        ax(col(6)).FontSize=6.8;
        %dscatter(airs_e,era5_a,'smoothing',2500,'bins',[50,50],'marker','o');
        colormap(ax(col(6)),CT)
        %ax(col(6)).CLim=[c_ll c_ul];
        clim([c_llimits(j) c_ulimits(j)])
        c=colorbar;
        c.FontSize=5.5;
        %c.Limits
        set(c,'TickDir','out');
        %if strcmp(field,'A')     
        if strcmp(field,'A')
            % c.TickLabels=num2cell([10^(-5) 10^(-4) 10^(-3)]);
            % c.Ticks=[10^(-5) 10^(-4) 10^(-3)];
            % c=colorbar(ax(col(6)),'YTick',log(c.Ticks),'YTickLabel',c.TickLabels);

            %  
            % %disp(c.Limits)
            c.Label.String='Fraction of log10(TC) x 10^{-3}';
            
            %c.TickLabels=num2cell(str2num(char(c.TickLabels))*(10^3));
            %c.Label.String='ND';
        else
            c.Ticks=[0 3e-3 6e-3 9e-3 12e-3];
            c.TickLabels=num2cell([0 3e-3 6e-3 9e-3 12e-3]);
            c.Label.String='Fraction of TC x 10^{-3}';
        end
        c.TickLabels=num2cell(str2num(char(c.TickLabels))*(10^3));
        c.Location='southoutside';
        axis square
        box off
        set(gca,'TickDir','out');
        set(gca,'Color',[0.88 0.88 0.88])
        %set(gca,'ColorScale','log')
        xlabel('AIRS')
        if j==1
            ylabel('ERA5 1')
        end
        xlim([xy_ll xy_ul])
        ylim([xy_ll xy_ul])
        label=l_col(6);
        st=subtitle(label);
        st.Units='normalized';
        st.Position=[0.115 0.83];
        hold on
        % 1:1 line
        plot(xy_ll:0.05:xy_ul,xy_ll:0.05:xy_ul,'Color',[0.5 0.5 0.5],'LineStyle','--')
        hold off

      
    end
    

    %t.Title.String=strjoin(day_night(i));

    exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\2D histograms\day_night plots 70th prc\' ...
                '2D_hist_', DN,'_log_10_cb.png'],''),'Resolution',600);
end 



