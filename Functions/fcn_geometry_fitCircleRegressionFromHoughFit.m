function [regression_fit_circle_center_and_radius, domain_box, radial_errors, standard_deviation]  =  ...
    fcn_geometry_fitCircleRegressionFromHoughFit(source_points, associated_points_in_domain, varargin)
% fcn_geometry_fitCircleRegressionFromHoughFit
% Given a set of points that are matched to a circle via a Hough vote,
% finds the circular regression fit circle and domain box. 
% 
% Format: 
% [regression_fit_circle_center_and_radius, domain_box] = fcn_geometry_fitCircleRegressionFromHoughFit(source_points,associated_points_in_domain, (fig_num))
%
% INPUTS:
%      source_points: a 3x2 matrix of the points used to create the Hough circle fit (used to find direction): 
% 
%      [point1_x  point1_y; point2_x  point2_y; point3_x  point3_y;]
%
%      These points are used to determine the direction of the resulting
%      circle fit.
%
%      associated_points_in_domain: an Nx2 list of points that should be
%      fit with regression, identified as within the domain according to
%      Hough voting.
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot the results.
%
% OUTPUTS:
%
%      regression_fit_circle_center_and_radius: a 3x1 matrix of the format:
%
%      [circleCenter_x  circleCenter_y; circleRadius]
%
%      domain_box: the box that encloses the 2-standard-deviation interval
%      around the regression circle fit.
%
%      radial_errors: the individual errors in each point, radially, in an
%      [N x 1] matrix
%
%      standard_deviation: the standard deviation in the errors
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%      fcn_geometry_fitSlopeInterceptNPoints 
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_fitCircleRegressionFromHoughFit
% for a full test suite.
%
% This function was written on 2024_01_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_09 - S Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
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
    debug_fig_num = 34838;
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
        narginchk(2,3);

        % Check the source_points input to be length exactly equal to 3
        fcn_DebugTools_checkInputsToFunctions(...
            source_points, '2column_of_numbers',[3 3]);

        % Check the associated_points_in_domain input to be length greater
        % than or equal to 3
        fcn_DebugTools_checkInputsToFunctions(...
            associated_points_in_domain, '2column_of_numbers',[3 4]);
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (3<= nargin) && (0==flag_max_speed)
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

% % Call a one-step circle fitting method
% [circleCenter, circleRadius] = fcn_INTERNAL_calc_Circle(associated_points_in_domain);

% Circle fit using Taubin's method
[circleCenter, circleRadius] = CircleFitByTaubin(associated_points_in_domain(:,1:2));

regression_fit_circle_center_and_radius = [circleCenter circleRadius];

% Find errors
radial_distances = sum((associated_points_in_domain-circleCenter).^2,2).^0.5;
radial_errors = radial_distances - circleRadius;
standard_deviation = std(radial_errors);

% Define the domain width
sigma_multiplier = 2;
max_orthogonal_distance = sigma_multiplier*standard_deviation; % max(abs(orthogonal_distances));

% Create a domain by doing a large range of angles across an inner and
% outer arc that spans the test area
angles = (0:1:360)'*pi/180;
inner_radius = max(0,(circleRadius - max_orthogonal_distance));
outer_radius = circleRadius + max_orthogonal_distance;
inner_arc = inner_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
outer_arc = outer_radius*[cos(angles) sin(angles)] + ones(length(angles(:,1)),1)*circleCenter;
domain_box = [inner_arc; flipud(outer_arc)];


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

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end        

    hold on;
    grid on;
    title('Regression fit of circle data');
    xlabel('X [meters]');
    ylabel('Y [meters]')
    
      
    % Plot the inputs: the source_points and associated_points_in_domain 
    plot(source_points(1,1),source_points(1,2),'g.','MarkerSize',30);
    plot(source_points(2,1),source_points(2,2),'b.','MarkerSize',30);
    plot(source_points(3,1),source_points(3,2),'r.','MarkerSize',30);
    h_plot = plot(associated_points_in_domain(:,1),associated_points_in_domain(:,2),'.','MarkerSize',10);
    current_color = get(h_plot,'Color');

    % Plot the circle fit
    fcn_geometry_plotCircle(circleCenter, circleRadius, current_color,fig_num);
    plot(circleCenter(1,1),circleCenter(1,2),'b+','MarkerSize',30);

    % Plot the domain
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
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

%% 
function [circleCenter, circleRadius] = fcn_INTERNAL_calc_Circle(associated_points_in_domain)
% Below uses method from:
% https://dtcenter.org/sites/default/files/community-code/met/docs/write-ups/circle_fit.pdf

mean_point = mean(associated_points_in_domain,1);
uv_points = associated_points_in_domain - mean_point;

Suu  = sum(uv_points(:,1).^2,1);
Suv  = sum(uv_points(:,1).*uv_points(:,2),1);
Svv  = sum(uv_points(:,2).^2,1);
Suuu = sum(uv_points(:,1).^3,1);
Svvv = sum(uv_points(:,2).^3,1);
Suuv = sum(uv_points(:,1).^2.*uv_points(:,2),1);
Suvv = sum(uv_points(:,1).*uv_points(:,2).^2,1);

A = [Suu Suv; Suv Svv];
b = 1/2*[Suuu+Suvv; Svvv+Suuv];
solution = A\b;

circleCenter = solution' + mean_point;
alpha = sum(solution.^2,1) + (Suu+Svv)/length(uv_points(:,1));
circleRadius = alpha^0.5;
end % Ends function


%%
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


%% 

% The following is wrong - not sure why?!
% The original (wrong) solution approach is based on the formula for a circle:
% 
% (x - x0)^2 + (y - y0)^2 = r0^2
% 
% where x0, y0 denotes the circle's center, and r0 is the radius. The goal
% is to find x0, y0, and r0 given data in the form of x,y.
%
% To solve this, we multiply out the squared terms for x,y, and rarrange
%
% (x^2 + y^2) - 2x*x0 - 2y*y0 + (x0^2+y0^2) = r0^2
%
% this is true for each x,y combination. Namely:
% 
% (x^2 + y^2)_1 - 2x_1*x0 - 2y_1*y0 + (x0^2+y0^2) = r0^2
%
% (x^2 + y^2)_2 - 2x_2*x0 - 2y_2*y0 + (x0^2+y0^2) = r0^2
%
% etc. where _1 and _2 (etc) represent the first point, the second point,
% etc.
% 
% If we take the difference between sequences of of the lines above, e.g.
% the data for (2-1), then:
%
% (x^2 + y^2)_2 - (x^2 + y^2)_1 - (2x_2 - 2x_1)*x0 - (2y_2 - 2y_1)*y0 = 0
% (x^2 + y^2)_3 - (x^2 + y^2)_2 - (2x_3 - 2x_2)*x0 - (2y_3 - 2y_2)*y0 = 0
%
% Rearranging, 
%
% (2x_2 - 2x_1)*x0 + (2y_2 - 2y_1)*y0 = ((x^2 + y^2)_2 - (x^2 + y^2)_1)
% (2x_3 - 2x_2)*x0 + (2y_3 - 2y_2)*y0 = ((x^2 + y^2)_3 - (x^2 + y^2)_2)
%
% which can be rewritten in linear regressor form as:
%
% A * X0 = b
%
% with:
%
% A = [
% (x_2 - x_1)  (y_2 - y_1);
% (x_3 - x_2)  (y_3 - y_2);
% (x_4 - x_3)  (y_4 - y_3);
% (etc)
% ];
%
% b = 1/2 * [
% ((x^2 + y^2)_2 - (x^2 + y^2)_1)
% ((x^2 + y^2)_3 - (x^2 + y^2)_2)
% ((x^2 + y^2)_4 - (x^2 + y^2)_3)
% (etc)
% ]
%
% and X0 = [x0; y0], or the circle center.
%
% The solution to this is given by the following steps:
%
%  AT*A * X0 = AT*b
%
%  X0 = inv(AT*A)*AT*b
%
%  or can be solved more efficiently by the pseudo inverse:
%
%  X0 = (AT*A)\AT*b
%   
% Once the center of the circle is known, the radius can be found by
% summing over all points the respective radial differences:
%
%  sum((x - x0)^2 + (y - y0)^2) = sum(r0^2)
%
% which gives:
%
% r0^2 = 1/N*sum((x - x0)^2 + (y - y0)^2)
% 
% 
% % Do some pre-calculations to fill in the A and b matricies
% radii_squared = sum(associated_points_in_domain.^2,2); % These are the radii-squared of points from origin, e.g. (x^2+y^2) for each point
% diff_points = diff(associated_points_in_domain,1,1);
% diff_rsquared = diff(radii_squared);
% 
% 
% % Construct the A-matrix and b matrix that will create the regressor.
% A = diff_points;
% b = 1/2*diff_rsquared;
% 
% 
% % Solve for the center points
% solution = (A'*A)\(A'*b);
% 
% % Reshape the centers into row format
% circleCenter = solution';
% 
% % NOTE: the following line is the slowest in the code. It can be sped
% % up if we do not take the square root
% radial_distances = sum((associated_points_in_domain-circleCenter).^2,2).^0.5;
% circleRadius = mean(radial_distances);
% radial_errors = radial_distances - circleRadius;
% standard_deviation = std(radial_errors);

