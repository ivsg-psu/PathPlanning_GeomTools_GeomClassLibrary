function domainShape = fcn_geometry_domainBoxByType(type_of_domain, varargin)   
% fcn_geometry_domainBoxByType -  generates the domain box given a
% particular type of domain. NOTE: the inputs change depending on the type
% of domain used.
%
% FORMAT:
%
%     domain_box = fcn_geometry_domainBoxByType(...
%     type_of_domain,...
%     (options),
%     (fig_num))
%
% INPUTS:
%
%      type_of_domain: a string indicating the type of domain to calculate.
%      Types include:
%
%          'arc': produces an arc type. If the inputs cause negative radii,
%          the arc is limited to zero radius. The call is:
%
%                domain_box = fcn_geometry_domainBoxByType(...
%                'arc',...
%                circleCenter, circleRadius, angles, distance_from_circle_to_boundary,
%                (fig_num))
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      domainShape: a polyshape created from the [x y] coordinates of the
%      bounding box representing the domain.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_domainBoxByType
% for a full test suite.
%
% This function was written on 2024_01_17 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2024_01_17 - S. Brennan
% -- wrote the function


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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

%% check input arguments?
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

if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(1,6); 
    
    % % Check the centers input
    % fcn_DebugTools_checkInputsToFunctions(...
    %     centers, '2column_of_numbers');
    %
    % % Use number of radii to calculate the number of centers
    % Ncircles = length(centers(:,1));
    %
    % % Check the radii input
    % fcn_DebugTools_checkInputsToFunctions(...
    %     radii, '1column_of_numbers',Ncircles);
    
end

% % Does user want to specify the degree_step?
% degree_step = 1;
% if 5 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         degree_step = temp;
%     end
% end
% 
% % Set plotting defaults
% plot_str = 'b-';
% plot_type = 1;  % Plot type refers to 1: a string is given or 2: a color is given - default is 1
% 
% % Check to see if user passed in a string or color style?
% if 6 <= nargin
%     input = varargin{2};
%     if ~isempty(input)
%         plot_str = input;
%         if isnumeric(plot_str)  % Numbers are a color style
%             plot_type = 2;
%         end
%     end
% end

% Does user want to specify the figure?
flag_do_plot = 1;
if (0==flag_max_speed) && (6==nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plot = 1;
    else
        flag_do_plot = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type_of_domain
    case {'arc'}
        circleCenter            = varargin{1};
        circleRadius            = varargin{2};
        angles                  = varargin{3};
        max_orthogonal_distance = varargin{4};

        inner_radius = max(0,(circleRadius - max_orthogonal_distance));
        outer_radius = circleRadius + max_orthogonal_distance;
        inner_arc = inner_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
        outer_arc = outer_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
        domain_box = [inner_arc; flipud(outer_arc)];
    otherwise
        error('Unknown domain type given: %s',type_of_domain)
end
domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);


%% Plot results?
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

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    axis equal
    grid minor;
    xlabel('X [meters]');
    ylabel('Y [meters]')

    % Plot the result
    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end