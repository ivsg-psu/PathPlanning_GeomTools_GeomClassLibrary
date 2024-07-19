function [curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = ...
fcn_geometry_curvatureAlongCurve(points_to_fit, varargin)
%% fcn_geometry_curvatureAlongCurve
% Given a set of XY data, finds the curvature along every point of the
% curve by finding the best "point" curvature at each XY location. The best
% curvature is defined as the variance-minimizing circle radius measured at
% that point, from points around that point. 
%
% The function calls curvatureAtPoint repeatedly, saving the best-fit
% curvatures, circle centers, signal to noise ratios, indicies for each
% best-fit region, etc. Because of the iterative nature of the function, it
% is very computationally heavy.
% 
% Format: 
% [curvatures, arc_centers, index_ranges, point_curvature_minimum, curvature_SNRs] = fcn_geometry_curvatureAlongCurve(points_to_fit, (data_width), (fig_num))
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
%      and right of the index_to_test. Default is to use all of them.
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      curvatures: the best-SNR curvature at each point. Note: curvature is
%      1/Radius.
%
%      arc_centers: the XY location at the best fit of the circle creating
%      the curvature at each point
%
%      index_ranges: the index range at each point that includes the
%      calculations for the curvature and SNR. The index range is the
%      number of indicies to left and right of the point that are included
%      in that point's curve fit
%
%      point_curvature_minimum: the smallest curvature that, with the
%      calculated tolerance of fit, would be indistinguishable from a line.
%      This is usefult to calculate the signal to noise ratio of the fit,
%      e.g. SNR = point_curvature/point_curvature_minimum
%
%      curvature_SNRs: the signal to noise ratio of the fit at each point,
%      equvalent to radius_maximum/arcRadius or
%      point_curvature/point_curvature_minimum. NOTE: these curvature SNRs
%      are corrected near the ends of the data to avoid incorrect
%      assessments, using data inset by the data_width if specified, or 20
%      points if not. If the true SNRs are needed, even at endpoints, they
%      can be calculated from the point_curvature_minimums.
%
% DEPENDENCIES:
%
%      fcn_geometry_curvatureAtPoint
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_curvatureAlongCurve
% for a full test suite.
%
% This function was written on 2024_07_18 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_07_18 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
        narginchk(1,3);

    end
end


% Does user want to specify data_width?
data_width = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        data_width = temp;
    end
end


% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (3<= nargin)
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
% How many points do we have?
Npoints = length(points_to_fit(:,1));

% Initialize storage arrays
index_curvatures = nan(Npoints,1); % Records which index was tested
curvatures = nan(Npoints,1);  % Records the best-fit curvature at each point
arc_centers = nan(Npoints,2); % Records the arc centers of each point
index_ranges = nan(Npoints,1); % Records how "wide" the curvature calculation was able to expand
point_curvature_minimums = nan(Npoints,1); % Records the largest radius (min curvature) that could possibly fit the curvature (best case)
curvature_SNRs = nan(Npoints,1); % Records the signal to noise ratio

% Initialize one of the plots
if flag_do_plots
    figure(fig_num);
    clf;

    subplot(1,3,1);
    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    plot(points_to_fit(:,1),points_to_fit(:,2),'k.','MarkerSize',20);

end

%% For each point, calculate curvature and SNR, and save results
for ith_point = 1:Npoints

    % Start by doing printing and drawing updates for this point
    if 1==flag_do_debug
        % Print results
        fprintf(1,'Testing point %.0d of %.0d\n',ith_point,Npoints);
    end
    % Draw position of test point
    if ith_point==1 && (1==flag_do_plots)
        figure(fig_num);
        h_testPoint = plot(points_to_fit(1,1),points_to_fit(1,2),'r.','MarkerSize',20);
    elseif (1==flag_do_plots)
        if Npoints<200 || (0==mod(ith_point,10))
            % Move to the current point
            set(h_testPoint,'XData', points_to_fit(ith_point,1),'YData',points_to_fit(ith_point,2));
            drawnow
        end
    end

    % Get the information for this test point
    [point_curvature, point_circle_center, index_range, point_curvature_minimum, best_SNR] = ...
        fcn_geometry_curvatureAtPoint(points_to_fit, ith_point, data_width, (-1));

    % Save results into the arrays
    index_curvatures(ith_point,1) = ith_point;
    curvatures(ith_point,1) = point_curvature;
    arc_centers(ith_point,:) = point_circle_center;
    index_ranges(ith_point,1) = index_range;
    point_curvature_minimums(ith_point,1) = point_curvature_minimum;
    curvature_SNRs(ith_point,1) = best_SNR;
end

%% Fix NaN values at start/end of data. 
% These are areas where the curvature
% is undefined since curves cannot be fit at the end of the data
% (regression requires minimum 3 points). To fix, we have the end points
% inheret the values from their closest neighbors.

N_to_fix = find(~isnan(curvatures),1)-1;

% Fix NaN values at start
index_to_inheret = N_to_fix+1;
range_to_fix = 1:N_to_fix;
curvatures(range_to_fix,1) = curvatures(index_to_inheret,1);
arc_centers(range_to_fix,1) = arc_centers(index_to_inheret,1);
arc_centers(range_to_fix,2) = arc_centers(index_to_inheret,2);
index_ranges(range_to_fix,1) = index_ranges(index_to_inheret,1);
point_curvature_minimums(range_to_fix,1) = point_curvature_minimums(index_to_inheret,1);
curvature_SNRs(range_to_fix,1) = curvature_SNRs(index_to_inheret,1);

% Fix NaN values at end
index_to_inheret = Npoints - N_to_fix;
range_to_fix = (Npoints - N_to_fix + 1):Npoints;
curvatures(range_to_fix,1) = curvatures(index_to_inheret,1);
arc_centers(range_to_fix,1) = arc_centers(index_to_inheret,1);
arc_centers(range_to_fix,2) = arc_centers(index_to_inheret,2);
index_ranges(range_to_fix,1) = index_ranges(index_to_inheret,1);
point_curvature_minimums(range_to_fix,1) = point_curvature_minimums(index_to_inheret,1);
curvature_SNRs(range_to_fix,1) = curvature_SNRs(index_to_inheret,1);

%% Fix SNR values at start/end of data. 
% The SNRs near the ends will be artifically low, and may result in bad
% data. To fix, we have the end points inheret the values from their
% closest neighbors.

if isempty(data_width)
    data_width = 20;
end

% Fix SNR values at start
index_to_inheret = data_width+1;
range_to_fix = 1:data_width;
curvature_SNRs(range_to_fix,1) = curvature_SNRs(index_to_inheret,1);

% Fix NaN values at end
index_to_inheret = Npoints - data_width;
range_to_fix = (Npoints - data_width + 1):Npoints;
curvature_SNRs(range_to_fix,1) = curvature_SNRs(index_to_inheret,1);


% Find the curvature SNR. First, we do not allow curavatures producing
% radii larger than 10000 meters
% point_curvature_minimums(curvatures<0.0001) = 1;
% curvature_SNRs = curvatures./point_curvature_minimums;





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
    % Plot the results

    %%%%%%%%%%%%%%%%%%%%%%
    figure(fig_num)
    subplot(1,3,2);
    cla;

    semilogy(index_curvatures,curvatures,'k-');
    hold on;
    semilogy(index_curvatures,point_curvature_minimums,'-','Color',[0.6 0.6 0.6]);

    grid on;
    xlabel('Index [count]');
    ylabel('Best curvature [1/m]');
    title('Curvatures')

    legend('Curvatures from points','Min curvature from noise')

    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,3,3);
    cla;

    semilogy(index_curvatures, curvature_SNRs,'k-');
    hold on;

    grid on;
    xlabel('Index [count]');
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

