function corrupted_points = fcn_geometry_corruptPointsWithOutliers(input_points, varargin)
% fcn_geometry_corruptPointsWithOutliers
% given a set of [Nx2] points, randomly adds outliers in the orthogonal
% direction using a random-normal magnitude.
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
%      outlier multiplied by a random-normal distribution. The default is
%      2.
%
%      fig_num: the figure number to use for plotting
%
% OUTPUTS:
%
%      corrupted_points: a list of test points that are corrupted
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_calcOrthogonalVector
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
% 2024_01_03 - S. Brennan
% -- added fast mode option
% -- added environmental variable options
% 2024_01_05 - S. Brennan
% -- switched to random-normal distributions
% -- fixed bug where only y values were being corrupted
% 2024_01_23 - S. Brennan
% -- added support for 3 dimensional points

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
        narginchk(1,4);

        % Check the points input to be length greater than or equal to 2
        fcn_DebugTools_checkInputsToFunctions(...
            input_points, '2or3column_of_numbers',[2 3]);

    end
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
if (0==flag_max_speed) && (4<= nargin)
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

    % Add random magnitudes onto orthogonal direction
    % y_range = max(input_points(:,2)) - min(input_points(:,2));
    % positive_or_negative = (rand(length(indicies_to_become_outliers),1)>0.5)*2.0 - 1;
    % magnitude_of_change = rand(length(indicies_to_become_outliers),1).*y_range.*magnitude_of_corruption./2 .* positive_or_negative;
    magnitude_of_change = randn(length(indicies_to_become_outliers),1).*magnitude_of_corruption;

    % Find unit vectors orthogonal to point-to-point vectors
    vectors_at_indicies = input_points(indicies_to_become_outliers,:) - input_points(indicies_to_become_outliers-1,:);
    orthogonal_unit_vectors_at_indicies = fcn_geometry_calcOrthogonalVector(vectors_at_indicies, [], -1); 
    
    % Corrupt using unit vectors
    corrupted_points(indicies_to_become_outliers,:) = corrupted_points(indicies_to_become_outliers,:) + magnitude_of_change.*orthogonal_unit_vectors_at_indicies;

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
    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end    

    hold on;
    grid on;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    if length(input_points(1,:))==2

        % Plot the input points
        plot(input_points(:,1),input_points(:,2),'k.','MarkerSize',20);

        % Plot the corrupted points
        plot(corrupted_points(:,1),corrupted_points(:,2),'m.','MarkerSize',15);
    else
        zlabel('Z [meters]')

        % Plot the input points
        plot3(input_points(:,1),input_points(:,2), input_points(:,3),'k.','MarkerSize',20);

        % Plot the corrupted points
        plot3(corrupted_points(:,1),corrupted_points(:,2),corrupted_points(:,3),'m.','MarkerSize',15);
    end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
    
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



