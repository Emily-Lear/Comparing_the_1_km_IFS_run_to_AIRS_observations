%% plot ERA5 winds for case studies

f_05='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\example case study ERA5 winds\05_ERA5_winds.mat';
f_09='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\example case study ERA5 winds\09_ERA5_winds.mat';

w_05=load(f_05);
w_09=load(f_09);


% find level at 39 km altitude 
alt= 39;
[~,level_i]=min(abs(w_05.z-alt));

u_05=w_05.u(:,:,level_i);
v_05=w_05.v(:,:,level_i);
u_09=w_09.u(:,:,level_i);
v_09=w_09.v(:,:,level_i);

%x,y]=ndgrid(w_05.x,w_05.y);
%% plot colour maps 

CT=flipud(cbrewer('div', 'RdBu', 30));

f=figure;

f.Position(1)=f.Position(1)-250;
f.Position(2)=f.Position(2)-400;
f.Position(3)=f.Position(3)*1; 
f.Position(4)=f.Position(4)*1.28; 

t=tiledlayout(2,2);
t.TileSpacing='tight';
%t.Padding='tight';
set(gcf,'color','w');

ax(1)=nexttile(1);
pcolor(w_05.x, w_05.y, u_05'); shading flat
axis equal
colormap(CT)
xlim([-1500 2000])
ylim([-2000 2000])
clim([-80 80])
set(gca,'TickDir','out');
set(gca,'Xticklabel',[])
t_l=text(0,0,'a)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];
title('U-component')

ax(2)=nexttile(2);
pcolor(w_05.x, w_05.y, v_05'); shading flat
axis equal
colormap(CT)
xlim([-1500 2000])
ylim([-2000 2000])
clim([-80 80])
set(gca,'TickDir','out');
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
t_l=text(0,0,'b)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];
title('V-component')

ax(3)=nexttile(3);
pcolor(w_09.x, w_09.y, u_09'); shading flat
axis equal 
colormap(CT)
xlim([-1500 2000])
ylim([-2000 2000])
clim([-80 80])
set(gca,'TickDir','out');
t_l=text(0,0,'c)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

ax(4)=nexttile(4);
pcolor(w_09.x, w_09.y, v_09'); shading flat
axis equal 
colormap(CT)
xlim([-1500 2000])
ylim([-2000 2000])
clim([-80 80])
set(gca,'TickDir','out');
c=colorbar;
c.Label.String= 'Wind speed (ms^{-1})';
c.Layout.Tile='east';
set(c,'TickDir','out');
set(gca,'Yticklabel',[])
t_l=text(0,0,'d)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

fontsize(10,"points")

% save figure
exportgraphics(gcf,'C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\example maps\winds and BT\winds_case_study.png','Resolution',300);



%% Plot wind contours 

f=figure;

f.Position(1)=f.Position(1)-250;
f.Position(2)=f.Position(2)-400;
f.Position(3)=f.Position(3)*1; 
f.Position(4)=f.Position(4)*1.28; 

t=tiledlayout(2,2);
t.TileSpacing='tight';
%t.Padding='tight';
set(gcf,'color','w');

ax(1)=nexttile(1);
[C,h]=contourf(w_05.x, w_05.y, u_05',[-20:5:80]); shading flat
set(h,'LineColor','none')
axis equal 
xlim([-1500 2000])
ylim([-2000 2000])
clim([-15 80])
set(gca,'TickDir','out');
set(gca,'Xticklabel',[])
t_l=text(0,0,'a)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

ax(2)=nexttile(2);
[C,h]=contourf(w_05.x, w_05.y, v_05',[-20:5:80]); shading flat
set(h,'LineColor','none')
axis equal
xlim([-1500 2000])
ylim([-2000 2000])
clim([-15 80])
set(gca,'TickDir','out');
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
t_l=text(0,0,'b)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

ax(3)=nexttile(3);
[C,h]=contourf(w_09.x, w_09.y, u_09',[-20:5:80]); shading flat
set(h,'LineColor','none')
axis equal 
xlim([-1500 2000])
ylim([-2000 2000])
clim([-15 80])
set(gca,'TickDir','out');
t_l=text(0,0,'c)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

ax(4)=nexttile(4);
[C,h]=contourf(w_09.x, w_09.y, v_09',[-20:5:80]); shading flat
set(h,'LineColor','none')
axis equal 
xlim([-1500 2000])
ylim([-2000 2000])
clim([-15 80])
set(gca,'TickDir','out');
c=colorbar;
c.Label.String= 'Wind speed (ms^{-1})';
c.Layout.Tile='east';
set(c,'TickDir','out');
set(gca,'Yticklabel',[])
t_l=text(0,0,'d)');
t_l.Units='normalized';
t_l.Position=[0.02 1.05];

fontsize(10,"points")



