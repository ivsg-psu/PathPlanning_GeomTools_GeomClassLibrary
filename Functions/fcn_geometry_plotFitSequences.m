function [XY_data, fit_numbers_of_data] = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters, varargin)
%% fcn_geometry_plotFitSequences
% Plots an individual geometry defined by a string name and parameter set.
%
% Format:
% XY_data = fcn_geometry_plotFitSequences(fitSequence_bestFitType, fitSequence_parameters, (segment_length), (format), (fig_num))
%
% INPUTS:
%
%      fitSequence_bestFitType: a cell array, in the order of fits, that
%      names the fit type.
% 
%      fitSequence_parameters: a cell array, in the order of fits, that
%      gives the parameters for each fit. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Vector regression segment fit'.
%
%      (OPTIONAL INPUTS)
%
%      segment_length: the smallest step to use for plotting, representing
%      the length (approximately) between points. Default is 0.1 meters.
%      This is a pass through parameter to fcn_geometry_plotGeometry - see
%      that function for examples.
%
%      format: A format string, e.g. 'b-', that dictates the plot style or
%      a color vector, e.g. [1 0 0.23], that dictates the line color. The
%      format string can also be a complex string - see the test script for
%      examples.
%      This is a pass through parameter to fcn_geometry_plotGeometry - see
%      that function for examples.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. If left empty, just plots to the current
%      figure.
%
% OUTPUTS:
%
%      XY_data: the data produced during plotting calculations. Note: this
%      data is returned even if fig_num is empty or set to -1.
%
%     fit_numbers_of_data: the specific fit, one value for each XY data
%     point, indicating which fit produced the XY data in the fit sequence
%
% DEPENDENCIES:
%      
%      fcn_geometry_plotGeometry
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_plotFitSequences
% for a full test suite.
%
% This function was written on 2024_04_14 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_04_14 - S. Brennan
% -- wrote the code
% 2024_06_26 - S. Brennan
% -- added segment_length and format input options 
% 2024_06_26 - S. Brennan
% -- added XY_data output
% 2024_07_20 - S. Brennan
% -- added model number output

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % Flag to plot the results for debugging
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
        narginchk(2,5);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '2column_of_numbers',[2 3]);
        %
        % % Check the transverse_tolerance input is a positive single number
        % fcn_DebugTools_checkInputsToFunctions(transverse_tolerance, 'positive_1column_of_numbers',1);
        %
        % % Check the station_tolerance input is a positive single number
        % if ~isempty(station_tolerance)
        %     fcn_DebugTools_checkInputsToFunctions(station_tolerance, 'positive_1column_of_numbers',1);
        % end
    end
end

% Does user want to specify the segment_length?
segment_length = [];
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        segment_length = temp;
    end
end

% Does user want to change the plot style?
% Set plotting defaults
plot_str = '';

% Check to see if user passed in a string or color style?
if 4 <= nargin
    input = varargin{2};
    if ~isempty(input)
        plot_str = input;
    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
if 0==flag_max_speed
    flag_do_plots = 1;
    if 5<= nargin
        temp = varargin{end};
        if ~isempty(temp)
            fig_num = temp;
        end
    else
        fig_num = gcf;
    end
end


%% Solve for the line fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_data = [];
fit_numbers_of_data = [];

% Fill in XY data from the domain fits
for ith_domain = 1:length(fitSequence_bestFitType)
    plot_type_string = fitSequence_bestFitType{ith_domain};
    parameters       = fitSequence_parameters{ith_domain};
    XY_domain = fcn_geometry_plotGeometry(plot_type_string, parameters, segment_length, plot_str, -1);
    XY_data   = [XY_data; XY_domain]; %#ok<AGROW>
    N_in_domain = length(XY_domain(:,1));
    fit_numbers_of_data = [fit_numbers_of_data; ith_domain*ones(N_in_domain,1)]; %#ok<AGROW>
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
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');

    
    % Plot the domain fits
    for ith_domain = 1:length(fitSequence_bestFitType)
        plot_type_string = fitSequence_bestFitType{ith_domain};
        parameters       = fitSequence_parameters{ith_domain};
        fcn_geometry_plotGeometry(plot_type_string, parameters, segment_length, plot_str, fig_num);
    end

    axis equal;
    
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


