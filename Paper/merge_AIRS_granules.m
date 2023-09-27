% Merge AIRS granules and model data resampled as AIRS into stipes

%% Append Granules into stripes of data % AIRS noise

% Directory 
D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise day_night';
DI_d=dir(D);
filenames_d=char(DI_d.name);
filenames_d(1:2,:)=[];

day=str2num(filenames_d(:,9:10));

% remove data after 10th Nov 2018
day_i=find(day>10);
filenames_d(day_i,:)=[];

% find difference between granule numbers of consequtive files 

% granule numbers
G_Num=str2num(filenames_d(:,19:21));

[Groups, Other]=consequtive_granules(G_Num);

% save list 
%save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\AIRS stripes.mat'), 'Groups');




%% loop through granules in stripes and save 

sz_g=size(Groups);


for j=1:sz_g(1)
    
    % structure for data
    variables=struct('time',[],'A',[],'HW',[],'VW',[],'MF',[],'AN',[] ...
        ,'Tp',[],'DN',[],'lon',[],'lat',[],'T',[],'MF_z',[],'MF_m',[]); % 
    granules=struct('day',[],'airs_time',[],'GN',[]);
    data=struct('ifs',variables,'s3d',variables,'era5',variables,'airs', ...
        variables,'x',[],'y',[],'point_spacing',[],'gn',granules);

    data_names=["ifs" "s3d" "era5" "airs"];
    v_names=["time" "A" "HW" "VW" "MF" "AN" "Tp" "DN" "lon" "lat" "T" "MF_z" "MF_m"]; %"MF_z" "MF_m"


    % loop through granules
    for file_i=Groups(j,1):Groups(j,2)

        filename=filenames_d(file_i,:);

        % load file
        f=fullfile(D,filenames_d(file_i,:));
    
      
        % add data to structure array
        g=load(f); 

        % if any points are within chosen region 
        if any(g.x(:)>=-4500 & g.x(:)<=4500 & g.y(:)>=-2000 & g.y(:)<=2000)

        data.gn.day=[data.gn.day string(filename(9:10))];
        data.gn.airs_time=[data.gn.airs_time string(filename(12:17))];
        data.gn.GN=[data.gn.GN string(filename(19:21))];

        data.x=[data.x g.x];
        data.y=[data.y g.y];
        data.point_spacing=g.ifs.point_spacing;

        % loop through datasets
        for m=1:length(data_names)
            name=data_names(m);

            % loop through variables
            for n=1:length(v_names)

            v_name=v_names(n);

            % append data from granule to stripe
            data.(name).(v_name)=[data.(name).(v_name) g.(name).(v_name)];          
            end

        end
        end

    end

    % if structure is not empty 
    if ~isempty(data.x)
        % save stripe (filename of last granule in stripe)
        save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night\data_ST_',data.gn.day(end), ...
            '_',data.gn.airs_time(end),'_',data.gn.GN(end),'.mat'),'-struct', 'data');
    end

end

%% load files not part of stripes and save to a file 
for i=1:length(Other)

     % file index
    file_i=Other(i);

    f=fullfile(D,filenames_d(file_i,:));
    data=load(f);
   
    % if any points are within chosen region 
    if any(data.x(:)>=-4500 & data.x(:)<=4500 & data.y(:)>=-2000 & data.y(:)<=2000)

    %save file 
    save(strcat('C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\granule stripes day_night\',filenames_d(file_i,:)),'-struct', 'data');

    end

end