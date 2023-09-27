% AIRS noise granules plot

% Directory 
d='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise granules\';

folders=["polar_night\" "midlat_night\" "trop_night\"];

% tile indicies 
TL=reshape(1:30,3,10);

f=figure;
t=tiledlayout(6,5); % 6 rows, 5 columns
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w');
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.63; %*0.4; %width %0.53
f.Position(4)=1170*0.7; 

% Colormaps
C=cbrewer('div', 'RdBu', 30);

% across-track and along track grid spacing 
x=0:19.95:89*19.95;
y=0:17.77:134*17.77;
[Y,X]=ndgrid(y,x);

for i=1:3

    folder=folders(i);
    d_folder=[strcat(d,folder)];
    DI_f=dir(d_folder);
    filenames=char(DI_f.name);
    filenames(1:2,:)=[];
    sz_f=size(filenames);

    % loop through granules
    for j=1:sz_f(1)

        filename=filenames(j,:);

        % filepath 
        filepath=fullfile(d_folder, filename);
        airs=load(filepath);

        ax(TL(i,j))=nexttile(TL(i,j));
        % plot airs granule at 39 km altitude 
        pcolor(X,Y,airs.Tp(:,:,5)'); shading flat 
        colormap(C)
        axis equal
        clim([-3 3])
        xlim([0 89*19.95])
        ylim([0 134*17.77])
        set(gca,'TickDir','out');
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        set(gca,'XTick',[]); set(gca,'YTick',[])

        Day=filename(15:16);
        GN=filename(23:25);
        Year=filename(20:21);

        title([Day '/11/' Year ' GN: ' GN ])

        if i==3 && j==10
            c=colorbar;
        end



    end
    

end 

title('AIRS Noise Granules')
fontsize(7,'points')

c.Label.String='Temperature Perturbation (K)';
set(c,'TickDir','out');
c.Layout.Tile='east';

%% Polar

% Directory 
d='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise granules\polar_night';
DI=dir(d);
filenames=char(DI.name);
filenames(1:2,:)=[];
sz_f=size(filenames);


f=figure;
t=tiledlayout(2,5); % 6 rows, 5 columns
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w');
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.8; %*0.4; %width %0.53
f.Position(4)=1170*0.3; 

% Colormaps
C=cbrewer('div', 'RdBu', 30);

% across-track and along track grid spacing 
x=0:19.95:89*19.95;
y=0:17.77:134*17.77;
[Y,X]=ndgrid(y,x);



% loop through granules
for j=1:sz_f(1)

    filename=filenames(j,:);

    % filepath 
    filepath=fullfile(d, filename);
    airs=load(filepath);

    ax(j)=nexttile(j);
    % plot airs granule at 39 km altitude 
    pcolor(X,Y,airs.Tp(:,:,5)'); shading flat 
    colormap(C)
    axis equal
    clim([-3 3])
    xlim([0 89*19.95])
    ylim([0 134*17.77])
    set(gca,'TickDir','out');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTick',[]); set(gca,'YTick',[])

    Day=filename(15:16);
    GN=filename(23:25);
    Year=filename(20:21);

    title([Day '/11/' Year ' GN: ' GN ])

    if j==10
       c=colorbar;
       set(c,'TickDir','out');
    end



end
    

title(t,'AIRS Noise Granules: Polar');

c.Label.String='Temperature Perturbation (K)';
c.Layout.Tile='east';
fontsize(9,'points')

% save figure
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\AIRS noise\AIRS noise regions\' ...
    'polar_AIRS_noise_granules.png'],'Resolution',300);

%% Midlatitudes

% Directory 
d='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise granules\midlat_night';
DI=dir(d);
filenames=char(DI.name);
filenames(1:2,:)=[];
sz_f=size(filenames);


f=figure;
t=tiledlayout(2,5); % 6 rows, 5 columns
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w');
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.8; %*0.4; %width %0.53
f.Position(4)=1170*0.3; 

% Colormaps
C=cbrewer('div', 'RdBu', 30);

% across-track and along track grid spacing 
x=0:19.95:89*19.95;
y=0:17.77:134*17.77;
[Y,X]=ndgrid(y,x);


% loop through granules
for j=1:sz_f(1)

    filename=filenames(j,:);

    % filepath 
    filepath=fullfile(d, filename);
    airs=load(filepath);

    ax(j)=nexttile(j);
    % plot airs granule at 39 km altitude 
    pcolor(X,Y,airs.Tp(:,:,5)'); shading flat 
    colormap(C)
    axis equal
    clim([-3 3])
    xlim([0 89*19.95])
    ylim([0 134*17.77])
    set(gca,'TickDir','out');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTick',[]); set(gca,'YTick',[])

    Day=filename(15:16);
    GN=filename(23:25);
    Year=filename(20:21);

    title([Day '/11/' Year ' GN: ' GN ])

    if  j==10
        c=colorbar;
    end



end
    

title(t,'AIRS Noise Granules: Midlatitudes')
fontsize(9,'points')

c.Label.String='Temperature Perturbation (K)';
set(c,'TickDir','out');
c.Layout.Tile='east';

% save figure
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\AIRS noise\AIRS noise regions\' ...
    'midlat_AIRS_noise_granules.png'],'Resolution',300);

%% Tropics

% Directory 
d='C:\Users\ejl45\OneDrive - University of Bath\AIRS\AIRS_noise\AIRS noise granules\trop_night';
DI=dir(d);
filenames=char(DI.name);
filenames(1:2,:)=[];
sz_f=size(filenames);


f=figure;
t=tiledlayout(2,5); % 6 rows, 5 columns
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w');
f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.8; %*0.4; %width %0.53
f.Position(4)=1170*0.3; 

% Colormaps
C=cbrewer('div', 'RdBu', 30);

% across-track and along track grid spacing 
x=0:19.95:89*19.95;
y=0:17.77:134*17.77;
[Y,X]=ndgrid(y,x);
   

% loop through granules
for j=1:sz_f(1)

    filename=filenames(j,:);

    % filepath 
    filepath=fullfile(d, filename);
    airs=load(filepath);

    ax(j)=nexttile(j);
    % plot airs granule at 39 km altitude 
    pcolor(X,Y,airs.Tp(:,:,5)'); shading flat 
    colormap(C)
    axis equal
    clim([-3 3])
    xlim([0 89*19.95])
    ylim([0 134*17.77])
    set(gca,'TickDir','out');
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTick',[]); set(gca,'YTick',[])

    Day=filename(15:16);
    GN=filename(23:25);
    Year=filename(20:21);

    title([Day '/11/' Year ' GN: ' GN ])

    if  j==10
        c=colorbar;
    end



end
    

title(t,'AIRS Noise Granules: Tropics')
fontsize(9,'points')

c.Label.String='Temperature Perturbation (K)';
set(c,'TickDir','out');
c.Layout.Tile='east';

% save figure
exportgraphics(gcf,['C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\AIRS noise\AIRS noise regions\' ...
    'trop_AIRS_noise_granules.png'],'Resolution',300);
