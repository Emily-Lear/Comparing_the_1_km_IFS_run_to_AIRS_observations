function DensityPlot=density_scatter(x,y,varargin)

% Function to plot scatter plots with density colour scale for points 

% INPUTS:
% x and y - data for scatterplot 

% OPTIONAL INPUTS:
% '2D_hist' - plot 2D histogram using histogram2 - 
% must be at the end of optional inputs.
% 'p_color'- plot the count density using pcolor instead of a scatter
% plot. - must be at the end of optional inputs.
% 'filled' - 'true' or 'false', whether the marker is filled or not.
% 'num_bins' - number of bins in x and y e.g. [50 50].
% 'mkr' - marker used for scatter plot e.g.'o'.
% 'sz' - size of points in scatter plot.
% 'xy_range' - set minimum and maximum bin edges on x and y axes 
% [min_x_value max_x_value min_y_value max_y_value].

% % Test:
% x=[1 2 3 4 5];
% y=[5 7 8 13 2];

% default values for optional inputs 
filled='true';
sz=10;
mkr='o';
num_bins=[50 50];
% % max and min values for bin edges
x_range=[min(x) max(x)];
y_range=[min(y) max(y)];

% loop through optional input names
for i=1:2:length(varargin)

    in_name=varargin{i};
    %in_val=varargin{i+1};


    switch(in_name)
        case 'filled'
            filled=varargin{i+1};
        case 'num_bins'
            num_bins=varargin{i+1};
        case 'mkr'
            mkr=varargin{i+1};
        case 'sz'
            sz=varargin{i+1};
        case 'xy_range'
            in_val=varargin{i+1};

            if length(in_val)==2
                x_range=in_val;
                y_range=in_val;
            elseif length(in_val)==4
                x_range=in_val(1:2);
                y_range=in_val(3:4);
            end
     end

end

% bin spacing 
x_spacing=(x_range(2)-x_range(1))./num_bins(1);
y_spacing=(y_range(2)-y_range(1))./num_bins(2);
% find bin edges
Xedges=linspace(x_range(1),x_range(2),num_bins(1)+1);
Yedges=linspace(y_range(1),y_range(2), num_bins(2)+1);


% find difference between min/max x/y values and min/max X/Y edges
% and find number of x/y spacings to include point
if min(x)<min(Xedges)
    n_min_x=ceil((min(Xedges)-min(x))/x_spacing);
else 
    n_min_x=1;
end
if max(x)>max(Xedges)
    n_max_x=ceil((max(x)-max(Xedges))/x_spacing);   
else 
    n_max_x=1;
end
if min(y)<min(Yedges)
    n_min_y=ceil((min(Yedges)-min(y))/y_spacing);
else 
    n_min_y=1;
end
if max(y)>max(Yedges)
    n_max_y=ceil((max(y)-max(Yedges))/y_spacing);   
else 
    n_max_y=1;
end


% new x and y range including all the points
x_range_new=[x_range(1)-(x_spacing*n_min_x) x_range(2)+(x_spacing*n_max_x)];
y_range_new=[y_range(1)-(y_spacing*n_min_y) y_range(2)+(y_spacing*n_max_y)];

% new x and y edges including all the data points 
Xedges=linspace(x_range_new(1),x_range_new(2),num_bins(1)+1);
Yedges=linspace(y_range_new(1),y_range_new(2), num_bins(2)+1);


% Xedges=[min(Xedges)-x_spacing Xedges max(Xedges)+x_spacing];
% Yedges=[min(Yedges)-y_spacing, Yedges, max(Yedges)+y_spacing];

% get bin counts 
[BinCounts, ~ ,~,binX, binY]=histcounts2(x,y,Xedges,Yedges);
ind=sub2ind(size(BinCounts),binX,binY);
% find density (bincounts/maximum)
CountDensity=BinCounts./max(BinCounts(:));

% find density value at each point in the scatterplot
c=CountDensity(ind);

if any(strcmpi(varargin,'2D_hist'))

    DensityPlot=histogram2(x,y,Xedges, Yedges,'DisplayStyle','Tile','Normalization','probability');
    % ,'Normalization','probability'
    shading flat;
    grid off
    box off

elseif any(strcmpi(varargin,'log_10_hist'))
    % log histogram count densities by 10
    %%% tick labels need to be changed back to actual values if this option
    %%% is used using formula:
    %%% log10(label_value)=plot_value-(|(min(log_CD))|+1)
    %%% Or: label_value=10^(plot_value-(|(min(log_CD))|+1))
    
    BC=BinCounts;
    BC(BC==0)=NaN;
    log_BC=log10(BC);
    log_BC(isnan(log_BC))=0;
    % can't plot negative values on the histogram so add 1+ (smallest log
    % value*-1)
    %CD_new=log_CD+abs(min(log_CD(:),[],'omitnan'))+1;
    DensityPlot=histogram2('XBinEdges',Xedges,'YBinEdges',Yedges,'BinCounts',log_BC,'DisplayStyle','Tile','Normalization','probability');

    disp(min(log_BC(:)))
    
    
elseif any(strcmpi(varargin,'p_color'))

    CountDensity(CountDensity==0)=NaN;        
    DensityPlot=pcolor(Xedges(2:end),Yedges(2:end),CountDensity'); shading flat;

elseif filled
% scatter plot

    DensityPlot=scatter(x,y,sz,c,mkr,'filled');
else 
    DensityPlot=scatter(x,y,sz,c,mkr);

end

% don't display output if nargout = 0
if nargout==0
    clear DensityPlot
end

end