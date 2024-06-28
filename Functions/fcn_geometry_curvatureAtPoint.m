function [point_curvature, point_circle_center, index_range, point_curvature_minimum] = fcn_geometry_curvatureAtPoint(points_to_fit, index_to_test, varargin)
%% fcn_geometry_curvatureAtPoint
% Given a set of XY data, finds the "point" curvature which is defined as
% the variance-minimizing circle radius measured at that point, from points
% around that point.
% 
% Format: 
% [point_curvature, radius_maximum] = fcn_geometry_curvatureAtPoint(points_to_fit, index_to_test, (fig_num))
%
% INPUTS:
%      points_to_fit: an [Nx2] matrix of N different [x y] points assumed
%      to be in sequence. Note: the function may break if the points are
%      not in sequence.
%
%      index_to_test: the index at which the curvature values should be
%      calculated.
%
%      (OPTIONAL INPUTS)
% 
%      data_width: the maximum number of index points to consider to left
%      and right of the index_to_test. Default is all of them.
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      point_curvature: the curvature of the minimum-variance circle fit.
%      note: curvature is 1/Radius.
%
%      point_circle_center: the XY location at the best fit of the circle
%
%      index_range: the index range of the points about the index_to_test
%      that includes the best-fit points used to calculate the circle fit
%
%      point_curvature_minimum: the smallest curvature that, with the
%      calculated tolerance of fit, would be indistinguishable from a line.
%      This is usefult to calculate the signal to noise ratio of the fit,
%      e.g. SNR = point_curvature/point_curvature_minimum
%
% DEPENDENCIES:
%
%      fcn_geometry_fitArcRegressionFromHoughFit
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_curvatureAtPoint
% for a full test suite.
%
% This function was written on 2024_06_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_06_27 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(2,4);

    end
end


% Does user want to specify data_width?
data_width = [];
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        data_width = temp;
    end
end

% % Does user want to specify flag_fit_backwards?
% flag_fit_backwards = 0;
% if (3<=nargin)
%     temp = varargin{2};
%     if ~isempty(temp)
%         flag_fit_backwards = temp;
%     end
% end



% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (4<= nargin)
    temp = varargin{end};
    if ~isempty(temp)        
        fig_num = temp;
        flag_do_plots = 1;
        
    end
end


%% Solve for the radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Fill in commonly used variables
Npoints = length(points_to_fit(:,1));

% Make sure that the index is not so close to the edge of the data that the
% regression will break. The regression needs 3 points, but there must be
% at least 5 points to calculate variance. The index point is centered
% between points, so the allowable indicies must be at least 2 indicies
% inward from each edge.

if index_to_test<3 || index_to_test>(Npoints-2)
    curvature_SNR = nan;
    index_width = nan;
    point_curvature = nan;
    point_circle_center = [nan nan];
    point_curvature_minimum = nan;
    index_range = nan;
else

    curvature_indicies = nan(Npoints,1);
    curvatures         = nan(Npoints,1);
    curvature_minimums = nan(Npoints,1);
    circle_centers     = nan(Npoints,2);
    if isempty(data_width)
        Nsearch = Npoints;
    else
        Nsearch = data_width;
    end
    
    for ith_expansion =2:Nsearch
        
        min_index = index_to_test-ith_expansion;
        max_index = index_to_test+ith_expansion;

        % Make sure the index range is inside range on indicies allowable
        % for points
        if min_index<=0
            break
        elseif max_index>Npoints
            break
        end
        [arcRadius, radius_maximum, arcCenter_xy] = fcn_geometry_fitRadiusToPoints(points_to_fit(min_index:max_index,:));
        curvature_indicies(ith_expansion,1) = ith_expansion;
        curvatures(ith_expansion,1)         = 1/arcRadius;
        curvature_minimums(ith_expansion,1) = 1/radius_maximum;
        circle_centers(ith_expansion,:)     = arcCenter_xy;

    end

    curvature_SNR = curvatures./curvature_minimums;
    [~,index_width] = max(curvature_SNR);

    point_curvature = curvatures(index_width);
    point_circle_center = circle_centers(index_width,:);
    point_curvature_minimum = curvature_minimums(index_width);
    index_range = index_width;
end





%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots

    
    %%%%%%%%%%%%%%%%%%%%%%
    figure(fig_num)
    subplot(1,3,1);
    cla;

    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    % Plot the input points
    plot(points_to_fit(:,1),points_to_fit(:,2),'b.','MarkerSize',20);

    % Make axis slightly larger?
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    temp_axis = axis;

    % Plot the circle fit at the point
    fcn_geometry_plotCircle(point_circle_center, 1/point_curvature,'r-',fig_num);
    
    
    % Plot the index range
    min_index = index_to_test-index_range;
    max_index = index_to_test+index_range;
    plot(points_to_fit(min_index:max_index,1),points_to_fit(min_index:max_index,2),'m.','MarkerSize',10)

    % Plot the query point
    plot(points_to_fit(index_to_test,1),points_to_fit(index_to_test,2),'g.','MarkerSize',30)

    axis(temp_axis);

    title('Input points');

    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,3,2);
    cla;

    semilogy(curvature_indicies,curvatures,'k-');
    hold on;
    semilogy(curvature_indicies,curvature_minimums,'-','Color',[0.6 0.6 0.6]);
    plot(index_width,point_curvature_minimum,'g.','Markersize',20);

    grid on;
    xlabel('index [count]');
    ylabel('curvature [1/m]');
    title('Curvatures')

    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,3,3);
    cla;
    grid on;
    hold on;

    plot(curvature_indicies, curvature_SNR,'k-');
    plot(curvature_indicies(index_width), curvature_SNR(index_width),'g.','Markersize',30);

    xlabel('index [count]');
    ylabel('SNR [unitless]');
    title('Curvature SNR')


end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function


%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

