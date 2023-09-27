function [Groups,Other] = consequtive_granules(Vector)

% Find groups of consequtive numbers in a list

% INPUT:
% Vector - vector of numbers

% OUTPUT:
% Groups -  array with index of 1st and last element of a group of
%           consequtive numbers in the input on each row
% Other -   array with index elements that are not part of the groups 
%           of consequtive numbers.

% % TEST: 
% % Directory 
% D='C:\Users\ejl45\OneDrive - University of Bath\ERA5, AIRS and IFS\39 km alt\wave properties new\data airs noise new 2';
% DI_d=dir(D);
% filenames_d=char(DI_d.name);
% filenames_d(1:2,:)=[];
% % granule numbers
% Vector=str2num(filenames_d(:,19:21));



GN_diff=diff(Vector); 

% if at end of day, change difference to 1 
GN_diff(GN_diff==-239)=1;

% find where diff is equal to 1 (index of consequtive number)
diff_1=find(GN_diff==1);
% index of 1s 
diff_2=diff_1+1;
ones_i=unique([diff_1 diff_2]);
% indicies of numbers in groups 
indicies=1:length(Vector);
indicies(ismember(indicies,ones_i))=[];
Other=indicies;

GN_1i=[];
% if difference is one and surrounding differences are not = 1
for i=2:length(GN_diff)
    if GN_diff(i-1)~=1 && GN_diff(i)==1 && GN_diff(i+1)~=1

        % get indicies 
        GN_1i=[GN_1i i];

    end
end

% end index of stripes 
GN_1end=GN_1i+1;



% loop through differences and remove values in middle of groups of 
% ones
for i=2:length(GN_diff)
    % if difference is one and surrounding differences are both ones 
    if GN_diff(i-1)==1 && GN_diff(i)==1 && GN_diff(i+1)==1
    
        % set value to nan
        GN_diff(i)=NaN;        
    end

    
end 

% get indicies of ones 
GN_i=find(GN_diff==1);

% indices of ones not surrounded by differences that are not ones
i2=find(~ismember(GN_i,GN_1i));
GN_i2=GN_i(i2);

%  end indices 
% reshape 
i_list=reshape(GN_i2,2,length(GN_i2)/2)';
% add 1 to second column to get index of last number in group
i_list(:,2)=i_list(:,2)+1;

% append 2 lists of start and end indices 
start_i=[GN_1i i_list(:,1)'];
[start_i, sI]=sort(start_i);
end_i=[GN_1end i_list(:,2)'];
end_i=end_i(sI);

% start and end values as 2 columns of an array
Groups=[start_i' end_i'];



end