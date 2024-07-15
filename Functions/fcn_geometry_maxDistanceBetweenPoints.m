function max_distance = fcn_geometry_maxDistanceBetweenPoints(points_to_check, varargin)
%% fcn_geometry_maxDistanceBetweenPoints
% Given a set of XY points, finds the maximum distance between the points,
% e.g. the maximum "span" of the points. The method is to find the
% Axis-Aligned Bounding Box (AABB) defining the maximum/minimum points in
% XY space, then find the diagonal. Note: this is not the TRUE maximum
% distance, but it will always upper-bound the actual maximum distance.
% 
% Format: 
% max_distance = fcn_geometry_maxDistanceBetweenPoints(points_to_check, (fig_num))
%
% INPUTS:
%      points_to_check: an [Nx2] matrix of N different [x y] points.
%
%      (OPTIONAL INPUTS)
% 
%      None
%
% OUTPUTS:
%
%      max_distance: the AABB-calculated maximum distance between points.
%
% DEPENDENCIES:
%
%      None
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_maxDistanceBetweenPoints
% for a full test suite.
%
% This function was written on 2024_06_27 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_06_27 - S. Brennan
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
        narginchk(1,2);

    end
end


% % Does user want to specify fitting_tolerance?
% fitting_tolerance = 0.1;
% if (2<=nargin)
%     temp = varargin{1};
%     if ~isempty(temp)
%         fitting_tolerance = temp;
%     end
% end
% 
% % Does user want to specify flag_fit_backwards?
% flag_fit_backwards = 0;
% if (3<=nargin)
%     temp = varargin{2};
%     if ~isempty(temp)
%         flag_fit_backwards = temp;
%     end
% end



% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if  (0==flag_max_speed) && (2<= nargin)
    temp = varargin{end};
    if ~isempty(temp)        
        fig_num = temp;
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
min_values = min(points_to_check,[],1);
max_values = max(points_to_check,[],1);
max_distance = real(sum((max_values-min_values).^2,2).^0.5);

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

    % Set up figure
    figure(fig_num);
    hold on;
    grid on;
    axis equal;
    xlabel('X [meters]');
    ylabel('Y [meters]');



    % Plot the AABB
    AABB_box = [...
        min_values(1,1) min_values(1,2);
        max_values(1,1) min_values(1,2);
        max_values(1,1) max_values(1,2);
        min_values(1,1) max_values(1,2);
        min_values(1,1) min_values(1,2);
        ];
    plot(AABB_box(:,1),AABB_box(:,2),'r-','LineWidth',3);

    % Plot the diagonal
    diagonal = [...
        min_values(1,1) min_values(1,2);
        max_values(1,1) max_values(1,2);
        ];
    plot(diagonal(:,1),diagonal(:,2),'m-','LineWidth',3);

    % Plot the input points very large
    current_color = fcn_geometry_fillColorFromNumberOrName(1,'points',[],-1);
    plot(points_to_check(:,1),points_to_check(:,2),'.','Color',current_color,'MarkerSize',10);

    % Put the title
    title(sprintf('Distance is: %.2f',max_distance));


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