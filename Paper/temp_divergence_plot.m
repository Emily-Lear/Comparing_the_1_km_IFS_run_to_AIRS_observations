% plot temperature divergence of the IFS resampled using both methods from 
% ERA5 and the IFS 


% scatter plot data noise
d='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\scatter region noise 1_14';

DI=dir(d);
filenames=char(DI.name);
filenames(1:2,:)=[];

day_prc=1:2:28;
night_prc=2:2:28;


% empty structure
stats=struct('corr',[],'P',[],'rmsd',[]);
day_night=struct('day',stats,'night',stats);
div=struct('ifs1_era5',day_night,'ifs2_era5',day_night,'ifs1_airs',day_night,...
    'ifs2_airs',day_night,'ifs1_ifs2',day_night,'era5_airs',day_night);

fields=["ifs1_era5" "ifs2_era5" "ifs1_airs" "ifs2_airs" "ifs1_ifs2" "era5_airs"];
names_1=["ifs" "s3d" "ifs" "s3d" "ifs" "era5"];
names_2=["era5" "era5" "airs" "airs" "s3d" "airs"];

for i=1:length(filenames(:,1))

    
    % find correlation and RMSD for each half day

    % load data
    f=fullfile(d,filenames(i,:));
    data=load(f);
        
    % daytime/ nighttime data
    if ismember(i,day_prc)
        DN='day';
    else
        DN='night';
    end

    

    for j=1:length(fields)


        field=fields(j);
        name_1=names_1(j);
        name_2=names_2(j);

        if i~=26
        
            % divergence of IFS from ERA5
            [C,P]=corrcoef(data.(name_1).T,data.(name_2).T,'Rows','complete');
            % add to list 
            div.(field).(DN).corr=[div.(field).(DN).corr C(2)];
            div.(field).(DN).P=[div.(field).(DN).P P(2)];
            % RMSD
            R=rmse(data.(name_1).T,data.(name_2).T,'omitnan');
            div.(field).(DN).rmsd=[div.(field).(DN).rmsd R]; 
            % don't include data for night 12 with missing AIRS data
        elseif i==26

        end
    end

end

% save data
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\temperature divergence\temp_div_1_14'),'-struct','div');

%% Plot temperature divergence
% 1st 14 days of November 2018

% Directory 
file='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\temperature divergence\temp_div_1_14.mat';
data=load(file);

day_night=["day" "night"];
dn=["Daytime" "Nighttime"];
fields=["ifs1_era5" "ifs2_era5" "ifs1_airs" "ifs2_airs" "ifs1_ifs2" "era5_airs"];
labels=["(a)" "(b)"];

% day / night 
for i=2 
    DN=day_night(i);

    f=figure;

    t=tiledlayout(1,2);
    t.TileSpacing='compact';
    %t.Padding='compact';
    set(gcf,'color','w')
        
    x=[1:12 14];
    % correlation
    ax(1)=nexttile(1);
    plot(x,data.ifs1_era5.(DN).corr,'-x','LineWidth',1)
    hold on 
    plot(x,data.ifs2_era5.(DN).corr,'-x','LineWidth',1)
    hold on
    plot(x,data.ifs1_airs.(DN).corr,'-x','LineWidth',1)
    hold on
    plot(x,data.ifs2_airs.(DN).corr,'-x','LineWidth',1)
    hold on
    plot(x,data.ifs1_ifs2.(DN).corr,'-x','LineWidth',1)
    hold on
    plot(x,data.era5_airs.(DN).corr,'-x','LineWidth',1)
    title('Correlation')
    ylabel('Correlation coefficient')
    xlim([1 14])
    ylim([0.7 1.02])
    t_l=text(0,0,labels(1));
    t_l.Units='normalized';
    t_l.Position=[0.015 1.07];
    t_l.FontSize=14;
    ax(1).FontSize=14;
    set(gca,'TickDir','out');
    hold off
    
    %rmsd
    ax(2)=nexttile(2);
    h1=plot(x,data.ifs1_era5.(DN).rmsd,'-x','LineWidth',1);
    hold on 
    h2=plot(x,data.ifs2_era5.(DN).rmsd,'-x','LineWidth',1);
    hold on
    h3=plot(x,data.ifs1_airs.(DN).rmsd,'-x','LineWidth',1);
    hold on
    h4=plot(x,data.ifs2_airs.(DN).rmsd,'-x','LineWidth',1);
    hold on
    h5=plot(x,data.ifs1_ifs2.(DN).rmsd,'-x','LineWidth',1);
    hold on
    h6=plot(x,data.era5_airs.(DN).rmsd,'-x','LineWidth',1);
    title('RMSD')
    ylabel('RMSD')
    xlim([1 14])
    ylim([0 8])
    t_l=text(0,0,labels(2));
    t_l.Units='normalized';
    t_l.Position=[0.015 1.07];
    t_l.FontSize=14;
    ax(2).FontSize=14;
    set(gca,'TickDir','out');
    hold off

    %title(t,dn(i))
    hLeg = copyobj([h1 h2 h3 h4 h5 h6],ax(2));
    set(hLeg, 'XData', NaN, 'YData', NaN, 'LineWidth', 1.5)
    leg=legend(hLeg,'IFS 1 & ERA5 1','IFS 2 & ERA5 1','IFS 1 & AIRS','IFS 2 & AIRS','IFS 1 & IFS 2','ERA5 1 & AIRS','Orientation','horizontal');
    leg.Layout.Tile='South';
    leg.FontSize=12;

 

    f.Position(1)=f.Position(1)-50;
    f.Position(2)=f.Position(2)-70;
    f.Position(3)=f.Position(3)*1.8;
    f.Position(4)=f.Position(4)*0.9;

    % save figure 
    exportgraphics(gcf,strjoin(['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\Temperature divergence\' ...
        'temp_div_1_14',DN,'_EGU.png'],''));

end







