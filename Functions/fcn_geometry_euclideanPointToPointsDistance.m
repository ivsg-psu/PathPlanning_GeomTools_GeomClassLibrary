function [dist] = ...
    fcn_geometry_euclideanPointToPointsDistance(...
    point1,...
    points2,...
    varargin)
% fcn_geometry_euclideanPointToPointsDistance calculates the 
% distance(s) between a single [1xd] point, POINT1, and another [Nxd] vector of
% points, POINTS2, where d is the dimension. Distance is returned as [Nx1]
% vector.
%
% FORMAT:
%
% [DIST] = fcn_geometry_euclideanPointToPointsDistance(POINT1, POINTS2, (fig_num))
%
% INPUTS:
%
%      POINT1: an 1x2 or 1x3 single xy or xyz point
%      in the form: [x1 y1 z1]
%
%      POINTS2: an Nx2 or Nx3 series of xy or xyz points 
%      in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      DIST: an N x  1 vector of distances [d1; d2; ... ; dn], where N is
%      the number of points in POINTS2, and DIST is the distance from each
%      of these points to POINT1
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%         pt1 = [1 1 5];
%         pt2 = [0 -3 -6; 34 1 17; 18 7 0];
%         dist=fcn_geometry_euclideanPointToPointsDistance(pt1,pt2);
%
% See the script: script_test_fcn_geometry_euclideanPointToPointsDistance
% for a full test suite.
%
% This function was originally written on 2018_11_17 by Seth Tau, supported
% now by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision History:
% 2021-05-28 - S. Brennan
% -- revised function to prep for geometry class 
% -- rewrote function to use vector sum
% -- added plotting option
% 2021-06-05
% -- fixed comments, added debugging option
% 2024_01_17 - Aneesh Batchu
% -- added max speed options 
% 2024_08_27 - Sean Brennan
% -- fixed function to clearly be different from PointsToPoints function

%% Debugging and Input checks

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

if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the point1 input
        fcn_DebugTools_checkInputsToFunctions(...
            point1, '2or3column_of_numbers', [1 1]);

        % Check the points2 input
        fcn_DebugTools_checkInputsToFunctions(...
            points2, '2or3column_of_numbers');
    end
end
    

% % Does user want to show the plots?
% if 3 == nargin
%     fig_num = varargin{end};
%     figure(fig_num);
%     flag_do_plot = 1;
% else
%     if flag_do_debug
%         fig = figure;
%         fig_num = fig.Number;
%         flag_do_plot = 1;
%     end
% end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (3 == nargin) 
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plot = 1;
    end
else
    if flag_do_debug
        fig = figure; 
        fig_num = fig.Number;
        flag_do_plot = 1;
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
Npoints = length(points2(:,1));

dist = sum((ones(Npoints,1)*point1-points2).^2,2).^0.5;

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

        
    midpoints = (point1+points2)/2;
    for ith_point=1:Npoints
        % 2D plot?
        if length(midpoints(1,:))==2
            % Plot the points
            xdata = [point1(1,1) points2(ith_point,1)];
            ydata = [point1(1,2) points2(ith_point,2)];
            plot(xdata,ydata,'.-','Linewidth',3,'Markersize',20);
            
            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),sprintf('d = %.1f',dist(ith_point,1)));
        else
            % Plot the points
            xdata = [point1(1,1) points2(ith_point,1)];
            ydata = [point1(1,2) points2(ith_point,2)];
            zdata = [point1(1,3) points2(ith_point,3)];
            plot3(xdata,ydata,zdata,'.-','Linewidth',3,'Markersize',20);

            % Label the midpoints
            text(midpoints(ith_point,1),midpoints(ith_point,2),midpoints(ith_point,3),sprintf('d = %.1f',dist(ith_point,1)));
            
            % Set to 3D view
            view(3);
        end
        
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

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function



