%% Plot topograpghy as regular distance

LonBox=[10 170];
LatBox=[0 90];
xRange=[-4700 4700];
yRange=[-2200 2200];
lon1=94;
lat1=52;
Loc1=[lat1 lon1];
point_spacing=[15 15 1 1];

[Topo_km] = topo_km_grid(LatBox,LonBox, xRange, yRange, Loc1, point_spacing);

% set O elevation to NaN
sea_i=find(Topo_km.elev <= 0);
Topo_km.elev2=Topo_km.elev;
Topo_km.elev2(sea_i)=NaN;


map1=[linspace(0,0,10)' linspace(0.25,0.35,10)' linspace(0.005,0.01,10)'];
map2=[linspace(0,1,120)' linspace(0.35,1,120)' linspace(0.01,1,120)'];
map3=[linspace(1,0.05,220)' linspace(1,0.05,220)' linspace(1,0.05,220)'];


map2(1,:)=[];
map3(1,:)=[];
map=[map1; map2; map3];

%Image 
figure
set(gcf,'color','w');

set(gca, 'visible', 'on')
pcolor(Topo_km.x, Topo_km.y, Topo_km.elev2')
set(gca,'Color','[0.1 0.55 0.8]')

shading flat
colormap(map);
axis equal
xlabel('Distance (km)')
ylabel('Distance (km)')
set(gca,'YDir','normal')
xlim([-4500 4500])
ylim([-2000 2000])
c=colorbar;
c.Label.String='Elevation (km)';
set(gca,'TickDir','out');
hold on
set(gca, 'TickDir', 'out')
contour(Topo_km.x, Topo_km.y, Topo_km.elev', [0,0],'-','color',[0 0 0]);


f=gcf;
exportgraphics(f,'file_path','Resolution',300)


% save structure 
save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\topography\Topo.mat'),'-struct', 'Topo_km');

        
