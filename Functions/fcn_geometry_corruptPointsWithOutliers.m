function corrupted_points = fcn_geometry_corruptPointsWithOutliers(input_points, varargin)
% fcn_geometry_corruptPointsWithOutliers
% given N points, with N>=2, creates a set of M points per unit distance
% between these points randomly distributed with variance sigma.
%
% corrupted_points = fcn_geometry_corruptPointsWithOutliers(input_points,
% (probability_of_corruption), (magnitude_of_corruption), (fig_num));
%
% INPUTS:
%
%      input_points: a Nx2 vector where N is the number of points, but at
%      least 2.
%
%      (Optional Inputs)
%
%      probability_of_corruption: the probabiity that a given point is an
%      outlier, from 0 to 1 (default is 0.02)
%
%      magnitude_of_corruption: the magnitude of corruption wherein the
%      outlier is sampled from a uniform distribution, as factor of the
%      y-axis range. The default is twice the y-range (2)
%
%      fig_num: the figure number to use for plotting
%
% OUTPUTS:
%
%      corrupted_points: a list of test points that are corrupted
%
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_corruptPointsWithOutliers
% for a full test suite.
%
% This function was written on 2023_12_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2023_12_12 
% -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(1,4);

    % Check the points input to be length greater than or equal to 2
    fcn_geometry_checkInputsToFunctions(...
        input_points, '2column_of_numbers',[2 3]);

end


% Does user want to specify probability_of_corruption?
probability_of_corruption = 0.02; % Default
if 2<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        probability_of_corruption = temp;
        if probability_of_corruption>1 || probability_of_corruption<0
            error('The probability_of_corruption must be between 0 and 1');
        end
    end
end

% Does user want to specify magnitude_of_corruption?
magnitude_of_corruption = 2; % Default
if 3<= nargin
    temp = varargin{2};
    if ~isempty(temp)
        magnitude_of_corruption = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if 4<= nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Corrupt the points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corrupted_points = input_points;

N_points = length(input_points(:,1));
random_flip = rand(N_points,1);
indicies_to_become_outliers = find(random_flip<=probability_of_corruption);

% Do not let first index be an outlier - it breaks the vector calculation
if ~isempty(indicies_to_become_outliers) 
    if indicies_to_become_outliers(1)==1
        indicies_to_become_outliers(1)=2;
    end

    % Calculate vectors in unit orthogonal direction
    vectors_at_indicies = input_points(indicies_to_become_outliers,:) - input_points(indicies_to_become_outliers-1,:);
    unit_vectors_at_indicies = fcn_geometry_calcUnitVector(vectors_at_indicies);
    orthogonal_unit_vectors_at_indicies = unit_vectors_at_indicies*[0 1; -1 0];

    % Add random magnitudes onto orthogonal direction
    y_range = max(input_points(:,2)) - min(input_points(:,2));
    positive_or_negative = (rand(length(indicies_to_become_outliers),1)>0.5)*2.0 - 1;
    try
        magnitude_of_change = rand(length(indicies_to_become_outliers),1).*y_range.*magnitude_of_corruption./2 .* positive_or_negative;
    catch
        disp('stop here')
    end

    corrupted_points(indicies_to_become_outliers,:) = corrupted_points(indicies_to_become_outliers,2) + magnitude_of_change.*orthogonal_unit_vectors_at_indicies;

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
    figure(fig_num);
    hold on;
    grid on;

    % Plot the input points
    plot(input_points(:,1),input_points(:,2),'g.','MarkerSize',20);
    
    % Plot the corrupted points
    plot(corrupted_points(:,1),corrupted_points(:,2),'m.','MarkerSize',15);

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



