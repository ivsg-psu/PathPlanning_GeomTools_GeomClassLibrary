function [point_density, counts1, counts2] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path, fig_histogram)
%% fcn_geometry_createHistogram
% Create a histogram that shows points per grid
% 
% FORMAT:
%
%      [point_density] = fcn_geometry_createHistogram(total_points_in_each_grid_in_the_driven_path, fig_histogram)
%
% INPUTS:     
%       
%      total_points_in_each_grid_in_the_driven_path: The total points in
%      each grid of the driven path of the mapping van
%
%      fig_histogram: the figure number of the histogram, (this is
%      mandatory input because the purpose of this function is to show
%      histogram plot).
%
% OUTPUTS:
%       
%      point_density: the density of point in each grid
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_geometry_createHistogram.m for a full
%       test suite.
%
% This code was written on 2024_07_18 by Aneesh Batchu
% This code was functionalized on 7/29/2024 by Jiabao Zhao
% Questions or comments? jpz5469@psu.edu
%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

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
%% Solve for the Maxs and Mins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Figure number of histogram
figure(fig_histogram); clf; 

% edges = (floor(min(total_points_in_each_grid_with_points_greater_than_zero)/10)*10):10:(ceil(max(total_points_in_each_grid_with_points_greater_than_zero)/10)*10); % Define the bin edges

% Create the histogram
% actual_driving_surface_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,edges,'Visible','on'); 
actual_driven_path_grids_hist = histogram(total_points_in_each_grid_in_the_driven_path,'Visible','on'); 
hold on 
% total_grids_hist = histogram(total_points_in_each_grid_with_points_greater_than_zero,edges,'Visible','on'); 
total_grids_greater_than_zero_hist = histogram(total_points_in_each_grid_in_the_driven_path,'Visible','on'); 

% Extract the counts for both histograms
counts1 = actual_driven_path_grids_hist.Values;
counts2 = total_grids_greater_than_zero_hist.Values;
binEdges = total_grids_greater_than_zero_hist.BinEdges;

% Calculate the overlapping counts
% overlapCounts = min(counts2, counts1);

% Find a ratio
point_density = sum(binEdges(1:2))/2; 

% Minimum number of points required 
% point_density =
% mean(total_points_in_each_grid_in_the_driven_path)/mean(total_points_in_each_grid_with_points_greater_than_zero);


% Add labels and title 
xlabel('Points per grid');
ylabel('Frequency');
title('Histogram of points per grid');
end
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
