function [parameters, standard_deviation_in_z, z_fit, unit_normal_vector, base_point, standard_deviation_in_plane_orthogonals, plane_distances, orthoBasisForPlane] = ...
    fcn_geometry_fitPlaneLinearRegression(points,varargin)
% fcn_geometry_fitPlaneLinearRegression
% Finds the coefficients, a, b, c... in the equation
%
%   ax + by + cz = Constant
%
% FORMAT: 
%
% [parameters, standard_deviation_in_z, z_fit, unit_normal_vector, base_point, standard_deviation_in_plane_orthogonals, orthoBasisForPlane] = ...
%    fcn_geometry_fitPlaneLinearRegression(points,(fig_num))
%
% INPUTS:
%
%      points: a Nx3 vector where N is the number of points, length N>=2.
%      Note that the plane fitting works in any dimension of 2 or higher
%
%      (OPTIONAL INPUTS)
% 
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      params: the [1 x 4] matrix of the parameters [a b c Constant]
%
%      standard_deviation_in_z: the standard deviation in the z-error of
%      the 
%  
%      z_fit: the model-fit z values
%
%      unit_normal_vector: the unit vector in the direction of <A, B, C> for the
%      equation: 
%
%             A*(x-x0) + B(y-y0) + C(z-z0) = 0
% 
%      base_point: the location closest to the origin of the plane
%
%      standard_deviation_in_plane_orthogonals: the standard deviation in
%      the point fitting error in the direction of the unit_normal_vector
%
%      plane_distances: the distances of each of the N points to the plane,
%      measured orthogonally from the plane. Returned as an [N x 1 ] vector.
%
%      orthoBasisForPlane: the orthonormal basis for the plane, e.g. the
%      "remainder" vectors that are within the plane created by excluding
%      the unit_normal_vector.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%       
% See the script: script_test_fcn_geometry_fitPlaneLinearRegression
% for a full test suite.
%
% This function was written on 2024_01_19 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2024_01_19 - S. Brennan
% -- wrote the code
% 2024_05_30 - S. Brennan
% -- rewrote the algorithm to allow vertical plane fits

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

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '3column_of_numbers',[3 4]);

    end
end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (2 == nargin) 
    temp_axis = varargin{1};
    if ~isempty(temp_axis)
        fig_num = temp_axis;
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
% What dimension is this?
dimension_of_points = length(points(1,:));

[normalToPlane, orthoBasisForPlane, base_point] = fcn_INTERNAL_affineFit(points);

magnitude_normal = sum(normalToPlane'.^2,2).^0.5;
unit_normal_vector = normalToPlane'./magnitude_normal;
Constant = sum(normalToPlane'.*base_point,2);
parameters = [normalToPlane' Constant];

z_planeFit = points(:,end);
N_points = length(points(:,1));

z_fit = (Constant*ones(N_points,1)  - points(:,1:end-1)*normalToPlane(1:end-1,1))./normalToPlane(end,1);
z_error = z_planeFit - z_fit;
standard_deviation_in_z = std(z_error, 0, 1, 'omitmissing');




%% Calculate the base point
% This is done by projecting all the input points via the unit vector
plane_distances = sum(ones(N_points,1)*unit_normal_vector.*points,2);
standard_deviation_in_plane_orthogonals = std(plane_distances);


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


    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end      
    
 
    % In subplot 1, plot the fit in 3D
    subplot(1,3,1);
    hold on;
    grid on;
    % axis equal;
    title('Plane fitting XYZ') 

    % Plot the points
    plot3(points(:,1),points(:,2),points(:,3),'k.','MarkerSize',20);
    plot3(points(:,1),points(:,2),z_fit,'.','MarkerSize',20);
    view(3);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.3;
        axis([temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,  temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y]);
    else
        temp_axis = axis;
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    ax = gca;
    ax.Clipping = "off";

    % Plot the plane
    reshaped_axis = reshape(temp_axis',2,[]);
    maxDifference = sum(diff(reshaped_axis).^2,2).^0.5;
    nudge = maxDifference*0.006;
    minX = reshaped_axis(1,1);
    maxX = reshaped_axis(2,1);
    minY = reshaped_axis(1,2);
    maxY = reshaped_axis(2,2);

    patchPoints_XY = [minX minY; maxX minY; maxX maxY; minX maxY];
    patchPoints_Z  = (Constant*ones(4,1)  - patchPoints_XY*normalToPlane(1:end-1,1))./normalToPlane(end,1);
    patchPoints = [patchPoints_XY patchPoints_Z];
    h_patch = patch(patchPoints(:,1), patchPoints(:,2), patchPoints(:,3), [0 0 1],'FaceAlpha',0.1); %#ok<NASGU>

    % x = [1 -1 -1 1]; % Generate data for x vertices
    x_planeFit = [temp_axis(2) temp_axis(1) temp_axis(1) temp_axis(2)]';

    % y = [1 1 -1 -1]; % Generate data for y vertices
    y_planeFit = [temp_axis(4) temp_axis(4) temp_axis(3) temp_axis(3)]';

    N_points = length(x_planeFit(:,1));
    z_planeFit = (Constant*ones(N_points,1)  - [x_planeFit y_planeFit]*normalToPlane(1:end-1,1))./normalToPlane(end,1); % Solve for z vertices data


    % Plot the base point
    plot3(base_point(1,1),base_point(1,2),base_point(1,3),'g.','MarkerSize',50);

    % Plot the unit vector
    quiver3(base_point(1,1),base_point(1,2),base_point(1,3), unit_normal_vector(1,1),unit_normal_vector(1,2),unit_normal_vector(1,3),0,'g','Linewidth',3);

    % Plot the vertical vector
    quiver3(x_planeFit(1,1),y_planeFit(1,1),z_planeFit(1,1), 0,0,1,0,'k','Linewidth',3);

    % Label the vertices with their numbers?
    if 1 == 1
        for ith_vertex = 1:length(points(:,1))
            if dimension_of_points==3
                text(points(ith_vertex,1)+nudge,points(ith_vertex,2),points(ith_vertex,3),...
                    sprintf('%.0d',ith_vertex),'Color',[0 0 0]);

            else
                text(points(ith_vertex,1)+nudge,points(ith_vertex,2),...
                    sprintf('%.0d',ith_vertex),'Color',[0 0 0]);
            end

        end
    end

    axis equal

    %% In subplot 2, plot the fit in 2D as X versus Z
    subplot(1,3,2);
    hold on;
    grid on;
    axis equal;
    title('Plane fitting X vs Z') 
    % Plot the points
    plot(points(:,1),points(:,3),'k.','MarkerSize',20);
    plot(points(:,1),z_fit,'.','MarkerSize',20);
   

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.3;
        axis([temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,  temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y]);
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');

    

    % Plot the "plane"
    plot(x_planeFit,z_planeFit,'-','Color',[0 0 1]);

    % Plot the base point
    plot(base_point(1,1),base_point(1,3),'g.','MarkerSize',50);

    % Plot the unit vector
    quiver(base_point(1,1),base_point(1,3), unit_normal_vector(1,1),unit_normal_vector(1,3),0,'g','Linewidth',3);

    % Plot the vertical vector
    quiver(x_planeFit(1,1),z_planeFit(1,1), 0,1,0,'k','Linewidth',3,'DisplayName','Vertical vector');

    axis equal

    %% In subplot 3, plot the fit in 2D as y versus z
    subplot(1,3,3);
    hold on;
    grid on;
    axis equal;
    title('Plane fitting Y vs Z') 
   % Plot the points
    plot(points(:,2),points(:,3),'k.','MarkerSize',20);
    plot(points(:,2),z_fit,'.','MarkerSize',20);
   

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.3;
        axis([temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,  temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y]);
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');


       % Plot the "plane"
    plot(y_planeFit,z_planeFit,'-','Color',[0 0 1]);

    % Plot the base point
    plot(base_point(1,2),base_point(1,3),'g.','MarkerSize',50);

    % Plot the unit vector
    quiver(base_point(1,2),base_point(1,3), unit_normal_vector(1,2),unit_normal_vector(1,3),0,'g','Linewidth',3);

    % Plot the vertical vector
    quiver(y_planeFit(1,1),z_planeFit(1,1), 0,1,0,'k','Linewidth',3);

    axis equal
    
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

%% OLD METHOD (WORKS, but not for vertical planes)
% 
% % The code below uses the least-y-squared solution of:
% % z = Ax + By + C
% %
% % z = [x y 1]*[A B C]'
% %
% % Let X = [x y 1], P = [A B C]'
% %
% % then
% %   z = X*P
% % 
% % which, solved for P:
% %
% % inv(X'*X)*(X'z) = P
% 
% x_planeFit = points(:,1);
% y_planeFit = points(:,2);
% N_points = length(x_planeFit(:,1));
% 
% X = [x_planeFit, y_planeFit, ones(length(points(:,1)),1)];
% z_planeFit = points(:,3);
% 
% if abs(det((X'*X)))>0.00001
%     P = (X'*X)\(X'*z_planeFit);
% else
%     P = nan(3,1);
% end
% 
% parameters = P;
% 
% z_fit = [x_planeFit y_planeFit ones(N_points,1)]*parameters; % Solve for z vertices data
% z_error = z_planeFit - z_fit;
% standard_deviation_in_z = std(z_error, 0, 1, 'omitmissing');
% 
% %% Calculate A, B, and C for unit vector
% C1 = parameters(1);
% C2 = parameters(2);
% 
% vector_projection = [-C1 -C2 1];
% 
% if ~any(isnan(vector_projection))
%     unit_normal_vector = fcn_geometry_calcUnitVector(vector_projection);
% else
%     unit_normal_vector = [nan nan nan];
% end
% 
% %% Calculate the base point
% % This is done by projecting all the input points via the unit vector
% plane_distances = sum(ones(N_points,1)*unit_normal_vector.*points,2);
% mean_plane_distance = mean(plane_distances);
% standard_deviation_in_plane_orthogonals = std(plane_distances);
% base_point = mean_plane_distance*unit_normal_vector;

%% fcn_INTERNAL_affineFit
% See: https://www.mathworks.com/matlabcentral/fileexchange/43305-plane-fit
function [n,V,p] = fcn_INTERNAL_affineFit(X)
%Computes the plane that fits best (lest square of the normal distance
%to the plane) a set of sample points.
%INPUTS:
%
% X: a N by 3 matrix where each line is a sample point
%
%OUTPUTS:
%
% n : a unit (column) vector normal to the plane
% V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
% plane
% p : a point belonging to the plane
%
% NB: this code actually works in any dimension (2,3,4,...)
% Author: Adrien Leygue
% Date: August 30 2013

%the mean of the samples belongs to the plane
p = mean(X,1);

%The samples are reduced:
R = bsxfun(@minus,X,p);
%Computation of the principal directions if the samples cloud
[V,~] = eig(R'*R);
%Extract the output from the eigenvectors
n = V(:,1);
V = V(:,2:end);

end % fcn_INTERNAL_affineFit