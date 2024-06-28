function [arcRadius, radius_maximum, circle_center_xy] = fcn_geometry_fitRadiusToPoints(points_to_fit, varargin)
%% fcn_geometry_fitRadiusToPoints
% Given a set of XY data, fits a circle to the data and finds the
% uncertainty in the fit, as a bounding box relative to 2 standard
% deviations.
% 
% Format: 
% [arcRadius, radius_maximum] = fcn_geometry_fitRadiusToPoints(points_to_fit, (fig_num))
%
% INPUTS:
%      points_to_fit: an [Nx2] matrix of N different [x y] points assumed
%      to be in sequence. Note: the function may break if the points are
%      not in sequence
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      arcRadius: the radius of the fit
%
%      radius_maximum: the maximum radius that would fit within the same
%      bounding box
%
%      circle_center_xy: the location of the best-fit circle center
%
% DEPENDENCIES:
%
%      fcn_geometry_fitArcRegressionFromHoughFit
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitRadiusToPoints
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
    debug_fig_num = 34838; 
else
    debug_fig_num = []; 
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


%% Solve for the radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% METHOD:
% 1: solve for the circle via regression fit
% 2: find the uncertainty in the fit
% 3: rotate the circle so that the box is AABB aligned
% 4: find the uncertainty in radius by box-fitting


% Fill in commonly used variables
Npoints = length(points_to_fit(:,1));

%% Step 0: Set up debugging figure?
if flag_do_debug
    figure(debug_fig_num);
    clf;
    subplot(3,2,1);
    cla;

    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    title('Input points');
    plot(points_to_fit(:,1),points_to_fit(:,2),'k.');

    % Make axis slightly larger?
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    temp_axis = axis;
end

%% Step 1: Circle fit using Taubin's method
[circle_center_xy, arcRadius] = CircleFitByTaubin(points_to_fit(:,1:2));

 

if flag_do_debug
    figure(debug_fig_num);
    subplot(3,2,2);
    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    title('Fit');
    plot(points_to_fit(:,1),points_to_fit(:,2),'k.');
    fcn_geometry_plotCircle(circle_center_xy, arcRadius,'r-',debug_fig_num);
    axis(temp_axis);
end


%% Step 2: Find standard deviation
% This is used later to find the bounding box
radial_distances = sum((points_to_fit-circle_center_xy).^2,2).^0.5;
radial_errors = radial_distances - arcRadius;
standard_deviation = std(radial_errors);

% Define the domain width
sigma_multiplier = 2;
max_orthogonal_distance = sigma_multiplier*standard_deviation; % max(abs(orthogonal_distances));

%% Step 3: Align the circle's box with middle of points
% Find the "arc middle" of points
vectors_from_center_to_points = points_to_fit - circle_center_xy;
unit_vectors_from_center_to_points = fcn_geometry_calcUnitVector(vectors_from_center_to_points);

% Which direction are the vectors oriented? Check with cross-product
cross_product_vectors = cross([unit_vectors_from_center_to_points(1:Npoints-1,:) zeros(Npoints-1,1)],[unit_vectors_from_center_to_points(2:Npoints,:) zeros(Npoints-1,1)],2);

% Remove any vectors that are repeated
vector_difference = sum(diff(unit_vectors_from_center_to_points).^2,2);
repeated_vector_indicies = vector_difference<1E-8;
cross_product_vectors(repeated_vector_indicies,:) = 0;

if all(cross_product_vectors(:,3)>=0) || all(cross_product_vectors(:,3)<=0)
    % Good data! 
else
    warning('on','backtrace');
    warning('The angles are not ordered. Unexpected results may occur as this function assumes only ordered points are used.')
end

% Is the arc clockwise or counter-clockwise?
flag_cross_product_is_positive = mean(cross_product_vectors(:,3))>=0;
if flag_cross_product_is_positive
    is_counterClockwise = 1;
    % Angles are all positive, e.g. the arc is CCW
    arc_angle_in_radians_start_to_end = fcn_geometry_findAngleUsing2PointsOnCircle(...
        circle_center_xy,...
        arcRadius,...
        points_to_fit(1,:),...
        points_to_fit(end,:),...
        is_counterClockwise, -1);
else
    is_counterClockwise = -1;
    % Angles are all positive, e.g. the arc is CCW
    arc_angle_in_radians_start_to_end = fcn_geometry_findAngleUsing2PointsOnCircle(...
        circle_center_xy,...
        arcRadius,...
        points_to_fit(1,:),...
        points_to_fit(end,:),...
        is_counterClockwise, -1);
end

half_angle = arc_angle_in_radians_start_to_end/2;
vector_center_to_start = points_to_fit(1,:)-circle_center_xy;
angle_start            = atan2(vector_center_to_start(2),vector_center_to_start(1));
if flag_cross_product_is_positive
    angle_middle           = angle_start + half_angle;
else
    angle_middle           = angle_start - half_angle;
end
position_middle = circle_center_xy + arcRadius*[cos(angle_middle) sin(angle_middle)];

if flag_do_debug
    
    figure(debug_fig_num);
    subplot(3,2,3);
    cla;
    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    title('Angle range');
    plot(points_to_fit(:,1),points_to_fit(:,2),'k.');
    points_start  = [circle_center_xy; points_to_fit(1,:)];
    points_middle = [circle_center_xy; position_middle];
    points_end    = [circle_center_xy; points_to_fit(end,:)];
    plot(points_start(:,1),points_start(:,2),'g-');
    plot(points_middle(:,1),points_middle(:,2),'b-');
    plot(points_end(:,1),points_end(:,2),'r-');

    axis(temp_axis);
end

%% Step 4: Find bounding box
box_width = arcRadius*abs(arc_angle_in_radians_start_to_end);
box_height = max_orthogonal_distance*2;
radius_maximum = fcn_geometry_maxRadiusInsideBox(box_width, box_height, -1);




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

    
    %%%%%%%%%%%%%%%%%%%%%%
    figure(fig_num)
    subplot(1,2,1);
    cla;

    hold on;
    grid on;
    axis equal
    xlabel('X [m]');
    ylabel('Y [m]');

    % Plot the input points
    plot(points_to_fit(:,1),points_to_fit(:,2),'k.');

    % Make axis slightly larger?
    temp = axis;
    %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

    temp_axis = axis;

    % Plot the circle fit
    fcn_geometry_plotCircle(circle_center_xy, arcRadius,'r-',fig_num);
    axis(temp_axis);

    % Plot angle range
    plot(points_to_fit(:,1),points_to_fit(:,2),'k.');
    points_start  = [circle_center_xy; points_to_fit(1,:)];
    points_middle = [circle_center_xy; position_middle];
    points_end    = [circle_center_xy; points_to_fit(end,:)];

    plot(points_start(:,1),points_start(:,2),'g-');
    plot(points_middle(:,1),points_middle(:,2),'b-');
    plot(points_end(:,1),points_end(:,2),'r-');

    %%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2);
    cla;

    hold on;
    grid on;
    axis equal
    xlabel('S [m]');
    ylabel('t [m]');

    % Plot the points in rotated format
    % Use the SE2 toolbox to transform
    % Start with translation of the arc's center to the origin

    translation_to_center        = -circle_center_xy;          % Push circle to be at origin
    transformMatrix_translation_into_center = se2(0,'theta',translation_to_center);

    % Rotate everything to push the start of arc angle to -90 -half_angle
    % degrees
    rotation_angle = -1*(angle_middle + pi/2);
    transformMatrix_rotation_into_St         = se2(rotation_angle,'theta',[0 0]);

    % Finally, translate everything upwards (counter-clockwise) or downwards (clockwise) to make arc1's end at the
    % origin
    translation_to_offset_center1  = [0 arcRadius];
    transformMatrix_offset_center1 = se2(0,'theta',translation_to_offset_center1);

    % Combine all the transformations
    St_transform  = transformMatrix_offset_center1*transformMatrix_rotation_into_St*transformMatrix_translation_into_center;

    % Move points
    points_rotated = transform(St_transform,points_to_fit);

    % Plot the bounding boxes
    fcn_geometry_maxRadiusInsideBox(2*abs(points_rotated(1,1)), abs(points_rotated(1,2)), fig_num);
        
    % Plot the input square
    square_points = [
        -1*box_width/2 0
        +1*box_width/2 0
        +1*box_width/2 box_height
        -1*box_width/2 box_height
        -1*box_width/2 0        
        ];
    plot(square_points(:,1),square_points(:,2),'g-','LineWidth',2);

    % Plot the rotated points
    plot(points_rotated(:,1),points_rotated(:,2),'r.');

    % Make axis slightly 
    temp = axis;
    axis_range_x = temp(2)-temp(1);
    axis_range_y = temp(4)-temp(3);
    percent_larger = 0.3;
    axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);




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

% The following was copied from fitCircleRegressionFromHoughFit

%% CircleFitByTaubin
% The following is Taubin's method of circle fitting, obtained from:
% https://www.mathworks.com/matlabcentral/fileexchange/22678-circle-fit-taubin-method?s_tid=prof_contriblnk
% This is a robust and accurate circle fit. It works well even if data
% points are observed only within a small arc. This circle fit was proposed
% by G. Taubin in article "Estimation Of Planar Curves, Surfaces And
% Nonplanar Space Curves Defined By Implicit Equations, With Applications
% To Edge And Range Image Segmentation", IEEE Trans. PAMI, Vol. 13, pages
% 1115-1138, (1991). It is more stable than the simple Circle Fit by Kasa
% (files
% #5557 and #22642) and slightly faster than Circle Fit by Pratt (file
% #22643).
function [circleCenter, circleRadius] = CircleFitByTaubin(XY)
%--------------------------------------------------------------------------
%  
%     Circle fit by Taubin
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With 
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------
n = size(XY,1);      % number of data points
centroid = mean(XY);   % the centroid of the data set
%     computing moments (note: all moments will be normed, i.e. divided by n)
Mxx = 0; Myy = 0; Mxy = 0; Mxz = 0; Myz = 0; Mzz = 0;
for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end

% Normalize
Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;


%    computing the coefficients of the characteristic polynomial
Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
A3 = 4*Mz;
A2 = -3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz - Mz*Mz*Mz;
A0 = Mxz*Mxz*Myy + Myz*Myz*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;
A33 = A3 + A3 + A3;
xnew = 0;
ynew = 1e+20;
epsilon = 1e-12;
IterMax = 20;
% Newton's method starting at x=0
for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
    if abs(ynew) > abs(yold)
       disp('Newton-Taubin goes wrong direction: |ynew| > |yold|');
       xnew = 0;
       break;
    end
    Dy = A1 + xnew*(A22 + xnew*A33);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Taubin will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Taubin negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end
%  computing the circle parameters
DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;
Par = [Center+centroid , sqrt(Center*Center'+Mz)];
circleCenter = Par(1,1:2);
circleRadius = Par(1,3);
end    %    CircleFitByTaubin
