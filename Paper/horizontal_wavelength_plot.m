% Horizontal wavelength - resampled models and AIRS 

% IFS 1 & 2

%Directory
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

% 70th percentile amplitude cutoffs 
f_prc='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\scatter plot data\day_night prc\day_night_prc_70.mat';
p=load(f_prc);

% % make stucture of indexes for each half day
file_indices=struct('d1',1,'d2',2:7, 'd3',8:14, 'd4',15:20, 'd5',21:26, ...
    'd6',27:32, 'd7',33:39, 'd8',40:46, 'd9',47:52, 'd10',53:58, ...
    'd11',59:64, 'd12',65:71, 'd13',72:77, 'd14',78:83, 'd15',84:89, ...
    'd16',90:96, 'd17',97:103, 'd18',104:109, 'd19',110:115, 'd20', ...
    116:121, 'd21',123:127);

half_days=fieldnames(file_indices);

% tiled layout indicies 
TL=reshape(1:40,4,10)';

% nighttime data for half days 
night=3:2:21;

% data names 
names=["airs" "ifs" "s3d" "era5"];

f=figure;
t=tiledlayout(10,4);
t.TileSpacing='tight';
t.Padding='compact';
set(gcf,'color','w');
t.Position=[0.1 0.13 0.85 0.83];

f.Position(2)=f.Position(2)-400;
f.Position(3)=830*0.7; %*0.4; %width %0.65
f.Position(4)=1170*0.63; % height 



i=0;

for d_i=night %11

    % tiled layout row
    i=i+1;
   
    % get file indexes 
    field=half_days{d_i};
    file_i=file_indices.(field);


    % data set
    for j=1:4
        
        %data name
        name=names(j);

        % loop through granules in night 
        for k=1:length(file_i)

            % file index
            index=file_i(k);
            % load file 
            f=fullfile(D,filenames_d(index,:));
            data=load(f);

            %smooth granule with 7x7 filter 
            s.(name).A=smoothn(data.(name).A,[7 7]);

            data.(name).HW(s.(name).A<p.(name).night)=NaN;

                     
            % save time for 1st granule
            if k==1
                if isfield(data,'gn')
    
                    day_1=char(data.gn.day(1));
                    time_1=char(data.gn.airs_time(1));            
                else                    
                     filename_1=filenames_d(index,:);
                     day_1=filename_1(9:10);
                     time_1=filename_1(12:17);  
                end
            end

            ax(TL(i,j))=nexttile(TL(i,j));
            pcolor(data.x, data.y, data.(name).HW); shading flat 
            % change colormap 
            CT=cbrewer('seq', 'Blues', 30);
            colormap(ax(TL(i,j)),CT);
            %c=colorbar;
            %c.Label.String='Horizontal Wavelength (km)';
            set(gca,'TickDir','out');
            set(gca,'Color',[0.8 0.8 0.8])
            clim([60 490]) 
            axis equal
            if i==1 && j==1
            title('AIRS')
            end
            if i==1 && j==2
                title('IFS 1')
            end
            if i==1 && j==3 
                title('IFS 2')
            end
            if i==1 && j==4 
                title('ERA5')
            end
            if i==10
            xlabel('Distance (km)')
            end
           
            
            xlim([-4500 4500])
            ylim([-2000 2000])
            hold(ax(TL(i,j)),'on') 


        end

        if i~=10
           set(gca,'XTickLabel',[]);
        end
        if j~=1
           set(gca,'YTickLabel',[]);
        end
        if j==1 
            ylabel('Distance (km)')
            st=subtitle(num2str(i),'FontWeight','bold');
            st.Units='normalized';
            st.Position;
            st.Position=[-0.45 0.37];
        % yl.Position(1)=yl.Position(1)-10;            
        end   

        

    end 


end

c=colorbar;
c.Label.String='Horizontal Wavelength (km)';
set(c,'TickDir','out');
c.Layout.Tile='south';

fontsize(7,'points')

% save figure
%exportgraphics(gcf,'C:\Users\ejl45\OneDrive - University of Bath\Plots\Nov 2018\ERA5_IFS_AIRS_colormaps\HW\nighttime_HW.png','Resolution',300);


