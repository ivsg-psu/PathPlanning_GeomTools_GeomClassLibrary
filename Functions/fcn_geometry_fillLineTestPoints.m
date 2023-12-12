
function test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, varargin)
% fcn_geometry_fillLineTestPoints
% given N points, with N>=2, creates a set of M points per unit distance
% between these points randomly distributed with variance sigma.
%
% [test_points] = fcn_geometry_fillLineTestPoints(seed_points, M, sigma)
%
% INPUTS:
%      seed_points: a Nx2 vector where N is the number of points, but at
%      least 2.
%
%      M: an integer listing the number of test points to generate per unit
%      distance.
%
%      sigma: athe standard deviation in points
%
% OUTPUTS:
%      test_points: a list of test points used to test regression fitting
%
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillLineTestPoints
% for a full test suite.
%
% NOTE: This function does NOT work for fitting all lines.
%
% This function was written on 2023_12_05 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_05 
% -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(3,4);

    % Check the points input to be length greater than or equal to 2
    fcn_geometry_checkInputsToFunctions(...
        seed_points, '2column_of_numbers',[2 3]);

end

% Does user want to show the plots?
if 4 == nargin
    fig_num = varargin{end};
    figure(fig_num);
    flag_do_debug = 1;
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
    end
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_segments = length(seed_points(:,1)) -1;

distances = sum((seed_points(2:end,:) - seed_points(1:end-1,:)).^2,2).^0.5;
unit_vectors = fcn_INTERNAL_calcUnitVector(seed_points(1:end-1,:),seed_points(2:end,:));
unit_orthogonals = unit_vectors*[0 1; -1 0];

test_points = [];
for ith_point = 1:length(distances)
    projection_distances = (0:(1/M):distances(ith_point))';
    N_points = length(projection_distances);
    orthogonal_distances = randn(N_points,1)*sigma;
    test_points = [test_points; seed_points(ith_point,:) + projection_distances*unit_vectors(ith_point,:) + orthogonal_distances*unit_orthogonals(ith_point,:)]; %#ok<AGROW>
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
if flag_do_debug
    figure(fig_num);
    hold on;
    grid on;

    % Plot the input points
    plot(seed_points(:,1),seed_points(:,2),'r.-','MarkerSize',20);
    
    % Plot the results
    plot(test_points(:,1),test_points(:,2),'b.','MarkerSize',10);

    % Make axis slightly larger
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    
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







%% fcn_INTERNAL_calcUnitVector
function unit_vector = fcn_INTERNAL_calcUnitVector(point_start,point_end)
vector_to_calculate    = point_end - point_start;
magnitude_vector_to_calculate = sum(vector_to_calculate.^2,2).^0.5;
unit_vector = vector_to_calculate./magnitude_vector_to_calculate;
end % Ends fcn_INTERNAL_calcUnitVector

