function [cone_parameters,fittedPoints, fitted_result] = fcn_geometry_fitRightCone(inputPoints,ring_id,varargin)

% fcn_geometry_fitRightCone calculates the ratio of radius to height of the
% cone
%
% FORMAT:
%
% cone_parameters = fcn_geometry_fitRightCone(inputPoints,ring_id)
%
% INPUTS:
%
%      inputPoints: a Nx3 vectors of point pairings in the one single scan
%      layer
%      
%      ring_id : an integer number indicates the id of the layer, ranging 
%      from 0 to 15
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      cone_parameters: an [1x2] vector containing ratio of radius to 
%       height of the cone and the id of the scan layer
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:   
%
% See the script: script_test_fcn_geometry_fitRightCone
% for a full test suite.
%
% This function was written on 2024_01_21 by X.Cao
% Questions or comments? xfc5113@psu.edu

% Revision history:
% 2024_01_21 - xfc5113@psu.edu
% -- original write of the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.

if (nargin==2)
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
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


if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the inputPoints input
        fcn_DebugTools_checkInputsToFunctions(...
            inputPoints, '3column_of_numbers');

       

end

% Does user want to show the plots?
flag_do_plots = 0;
if 3 == nargin
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

x_pts = inputPoints(:,1);
y_pts = inputPoints(:,2);
z_pts = inputPoints(:,3);
N_points = size(inputPoints,1);
R_square = x_pts.^2+y_pts.^2;
% calculates the c, the ratio of radius to height of the cone
c_ratio_square = mean((R_square)./(z_pts.^2));
c_ratio = sqrt(c_ratio_square);


% Estimate the Z with x, y and c
if ring_id<=7
    z_fitting = -sqrt(R_square./c_ratio_square);
else
    z_fitting = sqrt(R_square./c_ratio_square);
end
fittedPoints = [x_pts,y_pts,z_fitting];
cone_parameters = [c_ratio, ring_id];
fitted_error_sum = sum(vecnorm(inputPoints-fittedPoints,2,2));
fitted_error_mean = fitted_error_sum/N_points;
fitted_result = [fitted_error_sum, fitted_error_mean];
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
    axis equal;
    hold on;
    grid on;

   
    % Plot cone surface
    N_plot = 20;
    th_temp = linspace(0,2*pi,N_plot).';
    max_radius = 1.5*max(sqrt(inputPoints(:,1).^2+inputPoints(:,2).^2));
    
    radius_temp = linspace(0,max_radius,N_plot).';
    [R,Th] = meshgrid(radius_temp,th_temp) ;
    X_plot = R.*cos(Th) ;
    Y_plot = R.*sin(Th) ;
    if ring_id<=7
        Z_plot = -R/c_ratio;
    else
        Z_plot = R/c_ratio;
    end
    plot3(0,0,0,'m.','MarkerSize',20)
    hold on
    coneSurf = surf(X_plot,Y_plot,Z_plot);
    set(coneSurf,'FaceColor', [0 0 1], 'FaceAlpha',0.5)
     % Plot the input points
    plot3(inputPoints(:,1),inputPoints(:,2),inputPoints(:,3),'r.','MarkerSize',20);
    plot3(fittedPoints(:,1),fittedPoints(:,2),fittedPoints(:,3),'g.','MarkerSize',20)
    
    view(30,30)
    xlabel("X [m]",'fontsize',16)
    ylabel("Y [m]",'fontsize',16)
    zlabel("Z [m]",'fontsize',16)
    legend('Origin','Cone Surface','Input Points','Fitted Points','fontsize',16)
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

