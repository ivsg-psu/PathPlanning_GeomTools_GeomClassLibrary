function color_vector = fcn_geometry_fillColorFromNumberOrName(plot_number,varargin)
%% fcn_geometry_fillColorFromNumberOrName 
% Given a number, returns a color vector that for the chosen colorspace,
% always returns the same color. This is useful when there are many
% functions that are plotting the same data sequence onto different plots
% with different views or emphasis, and one wishes for all the plots to
% have the same color values. For example, all the plots with a plot_number
% of 2 will have the same color vector.
%
% The color space by default is 'gem12'.
% 
% Format: 
% color_vector = fcn_geometry_fillColorFromNumberOrName(plot_number, (string_identifier), (fig_num))
%
% INPUTS:
%      plot_number: an integer specifying the plot number to use for
%      retreiving a color. The colors appear in the sequence given by the
%      MATLAB command: "orderedcolors('gem12')". If the 'gem12' colorspace
%      is not available, then defaults to the "colororder" sequence, e.g.
%      the default sequence for MATLAB.
%
%      (OPTIONAL INPUTS)
% 
%      string_identifier: an optional string that is used to format plot
%      colors for specific common string types so that all have the same
%      color format. For example, the string types used for naming plotting
%      domain types for curve fits are supported, including the following:
%
%      'Regression arc': plots arcs in red
%
%      'Vector regression segment fit': plots lines or segments in blue
%      
%
%      fig_num: a figure number to plot results (not yet implemented). If
%      set to -1, skips any input checking or debugging, no figures will be
%      generated, and sets up code to maximize speed.
%
% OUTPUTS:
%
%      regression_domain: a structure that records details of the domain of
%      fitting. See fcn_geometry_fillEmptyDomainStructure for details.
%
%      std_dev_orthogonal_distance: the standard deviation in the point
%      fit, as measured in the transverse direction (orthogonal to the line
%      fit). E.g., this is the total-least-squares standard deviation.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_fitVectorToNPoints
%      fcn_geometry_fitSlopeInterceptNPoints 
%      fcn_geometry_domainBoxByType
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fillColorFromNumberOrName
% for a full test suite.
%
% This function was written on 2024_04_11 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_04_11 - S Brennan
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

        % % Check the source_points input to be length exactly equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     source_points, '2column_of_numbers',[2 2]);
        % 
        % % Check the associated_points_in_domain input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     associated_points_in_domain, '2column_of_numbers',[2 3]);
    end
end

% Does user want to specify best_fit_domain_box_projection_distance?
string_identifier = [];
if (2<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        string_identifier = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; %#ok<NASGU> % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (3<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp; %#ok<NASGU>
        flag_do_plots = 1;
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

% Get the color ordering?
try
    color_ordering = orderedcolors('gem12');
catch
    color_ordering = colororder;
end

N_colors = length(color_ordering(:,1));

if isempty(string_identifier)
    color_vector = color_ordering(mod(plot_number,N_colors)+1,:);
else
    switch lower(string_identifier)
        case {'arc','regression arc'}  % Arcs are red
            color_vector = [1 0 0];
        case {'line','segment','vector regression segment fit'} % Line fits are blue
            color_vector = [0 0 1];
        otherwise
            warning('on','backtrace');
            warning('Unrecognized plot color string: %s. Reverting to default.', string_identifier);
            color_vector = color_ordering(mod(plot_number,N_colors)+1,:);
    end
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
    % temp_h = figure(fig_num);
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end
    % 
    % % Get the color ordering?
    % try
    %     color_ordering = orderedcolors('gem12');
    % catch
    %     color_ordering = colororder;
    % end
    % 
    % N_colors = length(color_ordering(:,1));
    % 
    % hold on;
    % grid on;
    % axis equal;
    % 
    % % Plot the input points
    % plot(sorted_points_in_domain(:,1),sorted_points_in_domain(:,2),'k.','MarkerSize',10);
    % 
    % % Plot the fits    
    % ith_domain = 1;
    % color_vector = color_ordering(mod(ith_domain,N_colors)+1,:); %#ok<NASGU>
    % 
    % % Plot the base point and vector
    % plot(base_point_of_domain(1,1),base_point_of_domain(1,2),'g.','MarkerSize',30);
    % quiver(base_point_of_domain(:,1),base_point_of_domain(:,2),unit_tangent_vector_of_domain(:,1),unit_tangent_vector_of_domain(:,2),0,'g','Linewidth',5);
    % 
    % % Plot the domain
    % fcn_geometry_plotFitDomains(regression_domain, fig_num);
    % 
    % % Make axis slightly larger? And since this is the first one, save the
    % % axis limits.
    % if flag_rescale_axis
    %     temp = axis;
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     new_axis = [temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y];
    %     axis(new_axis);
    % end

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


