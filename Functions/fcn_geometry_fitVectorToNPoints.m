function [vector_root, unit_vector] = ...
    fcn_geometry_fitVectorToNPoints(points,varargin)
% fcn_geometry_fitVectorToNPoints
% Finds the vector fitting through a set of N points where N>=2. The
% solution performs a polar-form regression of the points, and thus works
% for point clusters in any orientation. It is based on the concept of
% averaging arc-tangent angles within the line segment after de-meaning the
% data. For more details, see:
%
%  "Feature Extraction and Scene Interpretation for Map-Based Navigation
%  and Map Building" by Kai Oliver Arras, Roland Y. Siegwart
%
% FORMAT: 
%
% [vector_root, unit_vector] = fcn_geometry_fitVectorToNPoints(points)
%
% INPUTS:
%
%      points: a Nx2 vector where N is the number of points, length N>=2. 
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      vector_root: the [2 x 1] matrix of the (x,y) point representing the
%      root of the vector, wherein the root is the point on the vector
%      closest to the origin.
%
%      unit_vector: the [2 x 1] matrix of the (deltax,deltay) length of the
%      unit vector attached to the root.
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%
% EXAMPLES:
%      
%      % BASIC example
%      points = [2 3; 4 5];
%      [vector_root, unit_vector] = fcn_geometry_fitVectorToNPoints(points)
% 
% See the script: script_test_fcn_geometry_fitVectorToNPoints
% for a full test suite.
%
% This function was written on 2021_05_26 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2020_06_25 
% -- wrote the code
% -- modified from fcn_geometry_fitVectorToNPoints
% 2024_01_13 - S. Brennan
% -- improved comments a bit
% -- added max speed options
% -- fixed bug wherein it crashes if empty figure given

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

        % Check the points input to be length greater than or equal to 2
        fcn_geometry_checkInputsToFunctions(...
            points, '2column_of_numbers',[2 3]);

    end
end

% Does user want to show the plots?
flag_do_plot = 0;
if (2 == nargin) && (0==flag_max_speed)
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
 
% The code below uses the polar form of the line, namely:
%
%     r_j = x cos(phi_j) + y sin(phi_j)
% 
% where the orthogonal distance of point (xi,yi) to this line is:
%     d_ortho(xi,yi) = r_j - xi*cos(phi_j) - yi*sin(phi_j)
%
% The regression solution is:
%
%  phi = 1/2*atan(-2sxy/(syy-sxx))
%
%  r = x_mean*cos(phi)+y_mean*sin(phi)
%
% with 
%    x_mean = mean(x), y_mean = mean(y), 
%    sxx = sum((xi - x_mean).^2)
%    syy = sum((yi - y_mean).^2)
%    sxy = sum((xi - x_mean)(yi - y_mean))

mean_point         = mean(points,1);
diff_mean          = points - mean_point.*ones(length(points(:,1)),1);
diff_squared_sum   = sum(diff_mean.^2,1);
diff_cross_sum     = sum(diff_mean(:,1).*diff_mean(:,2),1);

% Fill in variables from above 
x_bar = mean_point(1,1);
y_bar = mean_point(1,2);
sxx   = diff_squared_sum(1,1);
syy   = diff_squared_sum(1,2);
sxy   = diff_cross_sum;

% Find the solution for the polar form of the line
phi = 1/2*atan2(-2*sxy,syy-sxx);
r   = x_bar*cos(phi) + y_bar*sin(phi);

if r<0
    r = -1*r;
    phi = phi + pi;
end


% Convert polar form to vector form
vector_root = [r*cos(phi) r*sin(phi)];
unit_vector = [-sin(phi) cos(phi)]; % This is the perpendicular to the phi


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
    figure(fig_num);
    hold on;
    grid on;
    axis equal;
    
    plot(points(:,1),points(:,2),'r.');
    
    % Plot the root
    plot(vector_root(:,1),vector_root(:,2),'bo','Markersize',30);
    
    % Plot the unit vector
    quiver(vector_root(:,1),vector_root(:,2),unit_vector(:,1),unit_vector(:,2),0,'g','Linewidth',3);

    %     % Create an x vector of 100 points connecting lowest and highest points
    %     % in x, and plot the line fit
    %     x = linspace(min(points(:,1)),max(points(:,1)),100)';
    %
    %     % Fill in y avalues from y = (-A*x - C)/B;
    %     y = (x*(-A) -C*ones(length(x(:,1)),1))./B;
    %
    %     % If B is close to 0 (within tolerance of MATLAB), then the result is a
    %     % vertical line and we need to plot X = -C/A.
    %     if abs(B)<=eps(2)
    %         y = linspace(min(points(:,2)),max(points(:,2)),100)';
    %         x = -C/A.*ones(length(y),1);
    %     end
    %
    %     plot(x,y,'.-');
    
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


%% OLD METHOD (WRONG)

% % Find the mean values of the differences in the points, row-wise
% diff_points = diff(points);
% if length(diff_points(:,1))>1
%     mean_diff_points = mean(diff_points);
% else
%     mean_diff_points = diff_points;
% end
% 
% % Calculate x1, x2, y1, y2 by extracting each from points matrix
% x1 = points(1:end-1,1);
% y1 = points(1:end-1,2);
% x2 = points(2:end,1);
% y2 = points(2:end,2);
% 
% % Grab A, B, C values
% A = -1*mean_diff_points(1,2);
% B =    mean_diff_points(1,1);
% C = mean(x1.*y2 - y1.*x2);

