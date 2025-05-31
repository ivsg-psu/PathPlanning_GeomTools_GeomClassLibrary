function [unit_normal_vector, base_point, flags_in_agreement] = ...
    fcn_geometry_findPlaneNormal(points,varargin)
% fcn_geometry_findPlaneNormal
% Finds the the normal vector to a set of 3d points, keeping the
% cross-product direction. Note, this does not do regression fit (see
% fitPlaneLInearRegression), but rather takes the cross product of points
% in sequence and confirms that the resulting normal vectors agree with
% that of the first 3 points
%
% FORMAT: 
%
% [unit_normal_vector, base_point, flags_in_agreement] = ...
%    fcn_geometry_findPlaneNormal(points,(fig_num))
%
% INPUTS:
%
%      points: a Nx3 vector where N is the number of points, length N>=2.
%      Note that the plane fitting works in any dimension of 2 or higher
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      unit_normal_vector: the unit vector in the direction of <A, B, C> for the
%      equation: 
%
%             A*(x-x0) + B(y-y0) + C(z-z0) = 0
% 
%      base_point: the location centered on the points (mean of points)
%
%      flags_in_agreement: an Nx1 vector of 1 or 0, where N is the number
%      of points. A value of 1 is returned if the cross product of the
%      point agrees with the cross product of the first 3 points.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%       
% See the script: script_test_fcn_geometry_findPlaneNormal
% for a full test suite.
%
% This function was written on 2024_01_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2024_05_30 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        % Check the points input to be length greater than or equal to 3
        % rows, with 3 columns
        fcn_DebugTools_checkInputsToFunctions(...
            points, '3column_of_numbers',[3 4]);

    end
end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (2 == nargin) 
    temp_axis = varargin{1};
    if ~isempty(temp_axis)
        fig_num = temp_axis;
        figure(fig_num);
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number; %#ok<NASGU>
        flag_do_plot = 1;
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
% What dimension is this?
dimension_of_points = length(points(1,:));
assert(3==dimension_of_points);

base_point = mean(points,1);
vectors = diff(points,1,1);

unit_normal_vectors = cross(vectors(1:end-1,:),vectors(2:end,:),2);
unit_normal_vector = unit_normal_vectors(1,:);

flags_in_agreement = [1; 1; isequal(round(unit_normal_vector,6),round(unit_normal_vectors))];

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
if flag_do_plot


    % Print results
    fprintf(1,'\n\nResults of plane fitting: \n');
    table_data = [(1:length(points(:,1)))' [unit_normal_vector; unit_normal_vector; unit_normal_vectors] flags_in_agreement];
    header_strings = [{'Point No:'}, {'A'},{'B'},{'C'},{'Agree?'}]; % Headers for each column
    formatter_strings = [{'%.0f'},{'%.3f'},{'%.3f'},{'%.3f'},{'%.0f'}]; % How should each column be printed?
    N_chars = [12, 12, 12, 12, 12]; % Specify spaces for each column
    fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

    % Nothing to plot!
    
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
