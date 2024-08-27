function nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, varargin)
%fcn_geometry_anglesNearAngle  finds angles near a given angle
%
% Given a list of angles and an achor angle, finds the indicies of angles
% in the list that are within +/- angleRange of the anchorAngle and returns
% this as a [Mx1] vector where M is the number of indicies, sorted in the
% same order as anglesToSearch, that are within the angleRange. If no
% angles are within the range, returns an empty vector.
%
% The method to use is to convert all the anglesToSearch into unit vectors,
% and  also convert the anchorAngle to a unit vector. The dot product of
% these is then performed, and compared to cosine of the angleRange to find
% which vectors are "close". The reason for this method is to avoid
% wrap-around issues wherein angles such as 1 degree and 359 degrees are
% "close", or modulo issues wherein -179 degrees and 179 degrees are also
% close. The dot product method avoids both issues.
%
% FORMAT:
%
%       nearbyIndicies = fcn_geometry_anglesNearAngle(anchorAngle, anglesToSearch, angleRange, (fig_num))
%
% INPUTS:
%
%      anchorAngle: the angle, in radians, to search around
%
%      anglesToSearch: a [Nx1] vector of angles to check, in radians
%
%      angleRange: the range of angles around the anchorAngle to search, as
%      +/- range. For example, if an anchor angle is 2 radians and the
%      angleRange is 0.1 radians, then all angle indicies will be returned
%      that are between 1.9 and 2.1 radians.
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      nearbyIndicies: the indicies that have angles "nearby" the
%      anchorAngle
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%
%       script_test_fcn_geometry_anglesNearAngle
%
% This function was written on 2024_08_26 by Sean Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History
% 2024_08_26 S. Brennan
% -- started writing function

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS");
    MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG = getenv("MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PLOTCV2X_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PLOTCV2X_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
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

if 0 == flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(3,4);

    end
end

% Does user want to specify fig_num?
flag_do_plots = 0;
fig_num = []; % Initialize the figure number to be empty
if (0==flag_max_speed) && (4 <= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Write main code for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many angles are we searching
Nangles = length(anglesToSearch(:,1));

% Find unit vector for the anchorAngle
unitVector_anchorAngle= [cos(anchorAngle) sin(anchorAngle)];

% Find the unit vectors for all the anglesToSearch
unitVectors_anglesToSearch = [cos(anglesToSearch) sin(anglesToSearch)];


% Find dot products
dot_products = sum((ones(Nangles,1)*unitVector_anchorAngle).*unitVectors_anglesToSearch,2);

% Find angles in agreement - these are ones whose dot products
% are aligned, e.g. less than or equal to the cosine of the
% angle between the vectors
nearbyIndicies = find(dot_products>=cos(angleRange));


%% Any debugging?
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


if flag_do_plots == 1

    % Prepare the figure
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

    % Plot unit circle
    angles = (0:0.01:2*pi)';
    plot(cos(angles), sin(angles),'-','Color',[1 1 1]*0.5);

   % Plot the anchorAngle
   plot(unitVector_anchorAngle(:,1), unitVector_anchorAngle(:,2), 'b.','MarkerSize',20);

   % Plot the anglesToSearch
   plot(unitVectors_anglesToSearch(:,1),unitVectors_anglesToSearch(:,2),'k.','MarkerSize',5);

   % Plot the region
   anglesAlongRange = linspace(anchorAngle-angleRange,anchorAngle+angleRange,100)';
   regionBoundary = [0 0; cos(anglesAlongRange) sin(anglesAlongRange); 0 0];
   regionShape = polyshape(regionBoundary(:,1),regionBoundary(:,2),'Simplify',false,'KeepCollinearPoints',true);
   current_color = [0 1 0];
   plot(regionShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1,'EdgeAlpha',0);

   % Plot the points in agreement
   plot(unitVectors_anglesToSearch(nearbyIndicies,1), unitVectors_anglesToSearch(nearbyIndicies,2), 'go','MarkerSize',10,'LineWidth',2);


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

end
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
