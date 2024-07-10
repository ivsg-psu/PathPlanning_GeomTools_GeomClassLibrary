function fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, varargin)
%% fcn_geometry_printFitSequences
% Prints the parameters for a sequence of geometries to the console or a
% user-specified file ID.
%
% Format:
% fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))
%
% INPUTS:
%
%      fitSequence_fitTypes: a cell array, in the order of fits, that
%      names the fit types.
% 
%      fitSequence_parameters: a cell array, in the order of fits, that
%      gives the parameters for each fit. See
%      fcn_geometry_fillEmptyDomainStructure for details, specifically the
%      structure for 'Vector regression segment fit'.
%
%      (OPTIONAL INPUTS)
%
%      flag_all_parameters_same: sets the header/units printing style
%
%          if flag = 1, prints a descriptive header above all the
%          parameters, using the first parameter type as the format for all
%          of them. It then prints the parameters as numbers only (e.g.
%          2.4567). It prints the lead string and type ONLY on header line.
%          
% 
%          if flag = 0, then prints the parameters line by line,
%          as numbers then followed immediately by the units (for example:
%          2.4567m instead of 2.4567). Prints the lead string and type on
%          every line.
%
%          if flag is left empty, checks to see if all the parameters are
%          the same and sets flag=1 if they are, 0 otherwise. As well,
%          prints lead string and type on every line.
%
%      lead_string: A string that goes in front of the parameter set,
%      usually a descriptor.
%
%      fid: a file ID number to print results to a file. Default is 1 which
%      forces printing to the console.
%
% OUTPUTS:
%
%      XY_data: the data produced during plotting calculations. Note: this
%      data is returned even if fig_num is empty or set to -1.
%
% DEPENDENCIES:
%      
%      fcn_geometry_plotGeometry
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_printFitSequences
% for a full test suite.
%
% This function was written on 2024_07_10 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_07_10 - S. Brennan
% -- wrote the code using fcn_geometry_plotFitSequences as starter

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

% Does user want to specify the flag_all_parameters_same?
flag_all_parameters_same = [];
flag_forced_header = 0;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_all_parameters_same = temp;
        flag_forced_header = 1;
    end
end

if 0==flag_all_parameters_same
    flag_forced_header = 0;
end

% Does user want to specify the lead_string
lead_string = ' ';
if 4 <= nargin
    input = varargin{2};
    if ~isempty(input)
        lead_string = input;
    end
end

% Does user want to specify fid?
flag_do_plots = 0;
fid = 1; % Default is to print to console
if 0==flag_max_speed
    if 5<= nargin
        temp = varargin{end};
        if ~isempty(temp)
            fid = temp;
        end
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
Ndomains = length(fitSequence_fitTypes);

% If user did not specify print format, check to see if all the parameters
% are the same
if isempty(flag_all_parameters_same)
    flag_all_parameters_same = 1;
    for ith_domain = 1:Ndomains-1
        current_type_string = fitSequence_fitTypes{ith_domain};
        next_type_string = fitSequence_fitTypes{ith_domain+1};
        if ~strcmp(current_type_string,next_type_string)
            flag_all_parameters_same = 0;
        end
    end
end

% Print the header?
if 1==flag_all_parameters_same
    header_type_string = fitSequence_fitTypes{1};
    flag_print_header = 1;
    fcn_geometry_printGeometry(header_type_string, [], (flag_print_header), (lead_string), (fid))

end

% Continue printing lead string?
if 1==flag_forced_header
    lead_string = '';
end

% Print domain parameters

if 0==flag_forced_header
    flag_print_header = -1;
else
    flag_print_header = 0;
end

for ith_domain = 1:Ndomains
    header_type_string = fitSequence_fitTypes{ith_domain};
    fcn_geometry_printGeometry(header_type_string, fitSequence_parameters{ith_domain}, (flag_print_header), (lead_string), (fid))
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

% Nothing to plot
if flag_do_plots
    % temp_h = figure(fig_num);
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end
    % 
    % hold on;
    % grid on;
    % xlabel('X [m]');
    % ylabel('Y [m]');
    % 
    % 
    % % Plot the domain fits
    % for ith_domain = 1:length(fitSequence_fitTypes)
    %     plot_type_string = fitSequence_fitTypes{ith_domain};
    %     parameters       = fitSequence_parameters{ith_domain};
    %     fcn_geometry_plotGeometry(plot_type_string, parameters, segment_length, plot_str, fig_num);
    % end
    % 
    % axis equal;
    % 
    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
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


