function fcn_geometry_printGeometry(plot_type_string, parameters, varargin)
%% fcn_geometry_printGeometry
% Prints details of a geometry defined by a string name and parameter set.
% Default is to print results to the workspace, but optionally allows
% user-defined FID. 
%
% Format:
%      fcn_geometry_printGeometry(plot_type_string, parameters, (flag_print_header), (lead_string), (fig_num))
%
% INPUTS:
%
%      plot_type_string: a string indicating the geometry type to plot,
%      such as 'line', 'segment','spiral, or 'arc'. If plot string is
%      empty, or 'none', then nothing is plotted.
%
%      parameters: the parameter set describing the geometry. See
%      fcn_geometry_fillEmptyDomainStructure for details, as the parameter
%      set is different for each geometry type.
%
%      (OPTIONAL INPUTS)
%
%      flag_print_header: sets the header/units printing style
%
%          if flag=1, prints the descriptive header only without printing
%          the parameters (only need this at start of printing),
% 
%          if flag = 0 (default) then prints the parameters only as
%          numbers,
% 
%          if flag = -1 then prints the parameters as numbers then followed
%          immediately by the units (for example: 2.4567m instead of
%          2.4567).
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
%      (none)
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_printGeometry
% for a full test suite.
%
% This function was written on 2024_07_06 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2024_07_06 - S. Brennan
% -- wrote the code, using fcn_geometry_plotGeometry as a starter

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

% Does user want to specify the flag_print_header?
flag_print_header = 0;
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_print_header = temp;
    end
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

% Print the lead-in to the print
NleadCharacters = 20;

% Print the header string
fprintf(fid,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf(lead_string),NleadCharacters));

% Set up printing based on type
switch lower(plot_type_string)
    case {'none'}
        header_strings = {'(no parameters)'};
        simple_type = 'none';
        is_degrees  = 0;
        is_meters   = 0;

    case {'line'}
        header_strings = {'startX[m]','startY[m]','theta(deg)'};
        simple_type = 'line';
        is_degrees  = [0 0 1];
        is_meters   = [1 1 0];

    case {'segment','vector regression segment fit', 'line segment'}
        header_strings = {'startX[m]','startY[m]','theta(deg)','Slength[m]'};
        simple_type = 'segment';
        is_degrees  = [0 0 1 0];
        is_meters   = [1 1 0 1];

    case {'circle'}
        header_strings = {'centerX[m]','centerY[m]','radius[m]'};
        simple_type = 'circle';
        is_degrees  = [0 0 0];
        is_meters   = [1 1 1];

    case {'arc','regression arc'}
        header_strings = {'centerX[m]','centerY[m]','radius[m]','startAngle(deg)','endAngle(deg)','isCircle','CCW?'};
        simple_type = 'arc';
        is_degrees  = [0 0 0 1 1 0 0];
        is_meters   = [1 1 1 0 0 0 0];


    case {'spiral'}
        header_strings = {'startX[m]','startY[m]','startAngle(deg)','Slength(m)','K0(1/m)','Kf(1/m)'};
        simple_type = 'spiral';
        is_degrees  = [0 0 1 0  0  0];
        is_meters   = [1 1 0 1 -1 -1];

    otherwise
        warning('on','backtrace');
        warning('An error will now be thrown because a geometry string was not recognized.');
        error('Unknown plotting type: %s', plot_type_string);
end



% Print the header? If so, need to fill in the headers.
NumColumnChars = 10;

% Print the type
final_print_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s ',simple_type),NumColumnChars+1);


% Print header?
if 1==flag_print_header
    % Print headers. Note, the number of characters is +1 versus
    % specification because, for negative numbers (below), a space is
    % required before them and thus actual prints take +1 spaces.
    for ith_parameter = 1:length(header_strings)
        number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s ',header_strings{ith_parameter}),NumColumnChars+1);
        final_print_string = cat(2,final_print_string,number_string);

        % Do we need to print units? This adds 3 more spaces at end of each
        % entry
        if -1==flag_print_header
            final_print_string = cat(2,final_print_string,'   ');
        end
    end
    fprintf(fid,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',final_print_string),7*NumColumnChars));
else

    % Print parameters?
    if ~strcmp('none',simple_type)
        for ith_parameter = 1:length(header_strings)
            number_to_print = parameters(1,ith_parameter);

            % Do we need to print units? This adds three spaces to end of each
            % entry
            % Fill the default trailer string
            trailer_string = sprintf('%s','   ');
            if -1==flag_print_header
                if 1== is_meters(1,ith_parameter)
                    trailer_string = sprintf('%s','m  ');
                elseif -1== is_meters(1,ith_parameter)
                    trailer_string = sprintf('%s','1/m');
                elseif is_degrees(1,ith_parameter)
                    trailer_string = sprintf('%s','deg');
                end
            end

            % If the parameter is in degrees, make sure it is between 0 and 2pi
            % to start, then convert to degrees, and only print to 2 decimal
            % places. If it is not a degrees number, print to 4 decimal places.
            if 1==is_degrees(1,ith_parameter)
                number_to_print = mod(number_to_print,2*pi)*180/pi;
                number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.2f%s ',number_to_print,trailer_string),NumColumnChars);
            else
                number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.4f%s ',number_to_print,trailer_string),NumColumnChars);
            end

            

            % If start of the number has a minus sign, end-pad it. Otherwise,
            % front-pad it. This ensures the numbers have consistent numbers of
            % spaces.
            if strcmp(number_string(1),'-')
                final_print_string = cat(2,final_print_string,number_string, ' ');
            else
                final_print_string = cat(2,final_print_string,' ', number_string);
            end
        end
    end

    % The longest parameters ever printed require 7 columns
    fprintf(fid,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',final_print_string),7*NumColumnChars));
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
    % hold on;
    % grid on;
    % xlabel('X [m]');
    % ylabel('Y [m]');
    % 
    % % Get the color vector using the name
    % color_vector = fcn_geometry_fillColorFromNumberOrName(2,lower(plot_type_string));
    % 
    % % Check to see if need to ammend color specification to a complex plot_str
    % if length(plot_str)>3 && ~contains(lower(plot_str),'color')
    %     plot_str = cat(2,plot_str,',''Color'',color_vector');
    % end
    % 
    % % Plot lines as quiver arrows
    % if strcmp(plot_type_string,'line')
    %     direction_vector = XY_data(end,:)-XY_data(1,:);
    %     if plot_type==1
    %         if length(plot_str)>3
    %             eval_string = sprintf('quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2), 0, %s)',plot_str);
    %             eval(eval_string);
    %             eval_string = sprintf('quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2), 0, %s)',plot_str);
    %             eval(eval_string);
    %         elseif ~isempty(plot_str)
    %             quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,plot_str);
    %             quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,plot_str);
    %         else
    %             % Plot string is empty - use defaults
    %             quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,'-','LineWidth',3,'Color',color_vector);
    %             quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,'-','LineWidth',3,'Color',color_vector);
    %         end
    %     elseif plot_type==2
    %         quiver(XY_data(1,1),XY_data(1,2),    direction_vector(1,1),direction_vector(1,2),0,'-','LineWidth',3,'Color',plot_str);
    %         quiver(XY_data(end,1),XY_data(end,2),-direction_vector(1,1),-direction_vector(1,2),0,'-','LineWidth',3,'Color',plot_str);
    %     end
    % 
    % else
    %     % Plot everything else normally as XY data
    %     if plot_type==1
    %         if length(plot_str)>3                
    %             eval_string = sprintf('plot(XY_data(:,1),XY_data(:,2),%s)',plot_str);
    %             eval(eval_string);
    %         elseif ~isempty(plot_str)
    %             plot(XY_data(:,1),XY_data(:,2),plot_str);
    %         else
    %             % Plot string is empty - use defaults
    %             plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',color_vector);
    %         end
    %     elseif plot_type==2
    %         plot(XY_data(:,1),XY_data(:,2),'-','LineWidth',3,'Color',plot_str);
    %     end
    % 
    % end
    % 
    % % Plot green/red headers and tailers?
    % if ~any(isnan(XY_data),'all') && ~isempty(XY_data)
    %     if 1==1
    %         %%%%%%
    %         % Plot start and end locations as green/red points
    %         plot(XY_data(1,1),XY_data(1,2),     '.','Color',[0 1 0],'Linewidth',5,'MarkerSize',20);
    %         plot(XY_data(end,1),XY_data(end,2), 'o','Color',[1 0 0],'MarkerSize',10);
    % 
    %     else
    %         %%%%%%
    %         % Plot as start and end locations as green/red line segments
    % 
    %         % Plot green headers - calculated from vector direction
    %         maximum_arrow_length = 2*segment_length;
    %         minimum_arrow_length = 1*segment_length;
    % 
    %         vector_direction_start = XY_data(2,1:2) - XY_data(1,1:2);
    %         start_length = sum(vector_direction_start.^2,2).^0.5;
    %         unit_vector_direction_start = fcn_geometry_calcUnitVector(vector_direction_start);
    %         arrow_length = max(min(maximum_arrow_length,start_length*0.2),minimum_arrow_length);
    %         offset_start = XY_data(1,1:2) + arrow_length*unit_vector_direction_start;
    % 
    %         start_line = [XY_data(1,1:2) 0; offset_start, 0];
    %         plot(start_line(:,1),start_line(:,2), '-','Color',[0 1 0],'Linewidth',5);
    % 
    %         % Plot red tailers - calculated from vector direction
    %         vector_direction_end = (XY_data(end,1:2) - XY_data(end-1,1:2));
    %         end_length = sum(vector_direction_end.^2,2).^0.5;
    %         unit_vector_direction_end = fcn_geometry_calcUnitVector(vector_direction_end);
    %         arrow_length = max(min(maximum_arrow_length,end_length*0.2),minimum_arrow_length);
    %         offset_end = XY_data(end,1:2) - arrow_length*unit_vector_direction_end;
    %         end_line = [offset_end, 0; XY_data(end,1:2) 0];
    %         plot(end_line(:,1),end_line(:,2), '-','Color',[1 0 0],'Linewidth',5);
    %     end
    % end
    % 
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


