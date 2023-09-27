% Brightness temperature plots for example case studies 

f_nov_5="C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\example case study BT\airs_8um_bt_tibetan_plateau_20181105_for_emily.nc";

lon1=94;
lat1=52;

NCI=ncinfo(f_nov_5);
%NCI.Variables.Size;
Lon=double(ncread(f_nov_5,'lon'));
Lat=double(ncread(f_nov_5,'lat'));
BT=double(ncread(f_nov_5,'bt_8mu'));

% max(Lon(:))
% min(Lon(:))
% min(Lat(:))
% max(Lat(:))

lon=0:180;
lat=0:90;
[lat_1,lon_1]=ndgrid(lat,lon);
values=NaN(size(lon_1));

% regular distance grid 
x1=-1500:5:2000;
y1=-2000:5:2000;
[xq,yq]=ndgrid(x1,y1);

% find x and y coordinates from lon/lat values
Loc1=[lat1 lon1];
[x,y]=xy_distance(Lon,Lat,Loc1);


% interpolate to regular distance grid

F=scatteredInterpolant(x,y,BT,"nearest","none");
BT_grid=F(xq,yq);

% save BT data
data.BT=BT_grid;
data.x=xq;
data.y=yq; 

% save BT data
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\ERA5 BT example case study\ERA5_05_.mat'),'-struct', 'data');






%%
f_nov_9="C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\example case study BT\airs_8um_bt_tibetan_plateau_20181109_for_emily.nc";

lon1=94;
lat1=52;

Lon=double(ncread(f_nov_9,'lon'));
Lat=double(ncread(f_nov_9,'lat'));
BT=double(ncread(f_nov_9,'bt_8mu'));

% min(Lon(:))
% max(Lon(:))
% min(Lat(:))
% max(Lat(:))

lon=0:180;
lat=0:90;
[lat_1,lon_1]=ndgrid(lat,lon);
values=NaN(size(lon_1));

% regular distance grid 
x1=-1500:5:2000;
y1=-2000:5:2000;
[xq,yq]=ndgrid(x1,y1);

% find x and y coordinates from lon/lat values
Loc1=[lat1 lon1];
[x,y]=xy_distance(Lon,Lat,Loc1);


% interpolate to regular distance grid

F=scatteredInterpolant(x,y,BT,"nearest","none");
BT_grid=F(xq,yq);

data.BT=BT_grid;
data.x=xq;
data.y=yq; 

% save BT data
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\ERA5 BT example case study\ERA5_09.mat'),'-struct', 'data');

%% 5th & 9th Nov plot

f_05='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\ERA5 BT example case study\ERA5_05.mat';
f_09='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\ERA5 BT example case study\ERA5_09.mat';

d_05=load(f_05);
d_09=load(f_09);

BT2_05=d_05.BT; 
BT2_09=d_09.BT;
BT2_05(BT2_05>=220)=NaN;
BT2_09(BT2_09>=220)=NaN;

CT=cbrewer('div', 'Spectral', 30);

figure
t=tiledlayout(1,2);
t.TileSpacing='tight';
set(gcf,'color','w');

ax(1)=nexttile(1);
pcolor(d_05.x,d_05.y,d_05.BT); shading flat
colormap(flipud(CT))
axis equal
ylim([-2000 2000])
xlim([-1500 2000])
set(gca,'Color',[0.8 0.8 0.8])
set(gca,'TickDir','out');
clim([215 275])
hold on 
contour(d_05.x, d_05.y, d_05.BT, [220,220],'-','color',[1 0.4 0.8], 'LineWidth', 1);
hold off
t_l=text(0,0,"a)");
t_l.Units='normalized';
t_l.Position=[0.02 1.07];
t_l.FontSize=10;

ax(2)=nexttile(2);
pcolor(d_09.x,d_09.y,d_09.BT); shading flat
colormap(flipud(CT))
axis equal
ylim([-2000 2000])
xlim([-1500 2000])
set(gca,'Color',[0.8 0.8 0.8])
set(gca,'TickDir','out')
clim([215 275])
set(gca,'Yticklabel',[])
hold on 
contour(d_09.x, d_09.y, d_09.BT, [220,220],'-','color',[1 0.4 0.8], 'LineWidth', 1);
hold off
c=colorbar;
c.Label.String= '8 μm Brightness temperature (K)';
set(c,'TickDir','out');
t_l=text(0,0,"b)");
t_l.Units='normalized';
t_l.Position=[0.02 1.07];
t_l.FontSize=10;

% save plot 
exportgraphics(gcf,'C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\example maps\BT\BT_case_study.png','Resolution',300);



