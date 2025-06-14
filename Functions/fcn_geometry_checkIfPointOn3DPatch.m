function [flag_isInside3Dpatch, flag_isOn3DPatchPlane, flag_projectsInsidePatch] = ...
    fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, varargin)
% fcn_geometry_checkIfPointOn3DPatch
% Finds the coefficients, a, b, c... in the equation
%
%   ax + by + cz = Constant
%
% FORMAT: 
%
% [flag_isInside3Dpatch, flag_isOn3DPatchPlane] = fcn_geometry_checkIfPointOn3DPatch(patchPoints, testPoints, (tolerance), (fig_num))
%
% INPUTS:
%
%      patchPoints: a Nx3 vector where N is the number of points, length N>=3
%      that define which points define the patch. Note: the patch can be
%      nonconvex within its plane, but all points must be on the same
%      plane.
%
%      testPoints: a Mx3 vector where M is the number of testing points,
%      length M>=1 that define which points are to be tested to see if they
%      are both in the patchPoint plane and also bounded by these points.
%
%      (OPTIONAL INPUTS)
% 
%      tolerance: the tolerance to test if points are within the plane and
%      near edges. Default is eps*1000;
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      flag_isInside3Dpatch: an Mx1 array of 1's or 0's, with 1 indicating
%      that the respective test point is within the patchPoint plane and
%      enclosed by the patchPoint boundary. Note: flag_isInside3Dpatch is
%      simply the element-wise multiplication of flag_isOn3DPatchPlane and
%      flag_projectsInsidePatch. 
%
%      flag_isOn3DPatchPlane: an Mx1 array of 1's or 0's, with 1 indicating
%      that the respective test point is within the patchPoint plane.
%
%      flag_projectsInsidePatch: an Mx1 array of 1's or 0's, with 1 indicating
%      that the respective test point is enclosed by the 3D patch, when
%      projecting that point onto the plane created by the 3D patch.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%       
% See the script: script_test_fcn_geometry_checkIfPointOn3DPatch
% for a full test suite.
%
% This function was written on 2025_06_12 by S. Brennan
% Questions or comments? sbrennan@psu.edu 


% Revision history:
% 2025_06_12 - S. Brennan
% -- wrote the code


%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(2,4);

        % % Check the points input to be length greater than or equal to 2
        % fcn_DebugTools_checkInputsToFunctions(...
        %     points, '3column_of_numbers',[3 4]);

    end
end

% Does user want to specify tolerance?
tolerance = 1000*eps;
if (3<=nargin)
    temp = varargin{1};
    if ~isempty(temp)
        tolerance = temp;
    end
end

% Does user want to show the plots?
flag_do_plot = 0;
if (0==flag_max_speed) && (4 == nargin) 
    temp = varargin{end};
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

% Make sure points are in 3D!
dimension_of_points = length(patchPoints(1,:));
if dimension_of_points~=3
    warning('on','backtrace');
    warning('3d points were expected, but not given');
    error('Expected 3D points for plane fitting');
end

if dimension_of_points~=length(testPoints(1,:))
    warning('on','backtrace');
    warning('3d points were expected, but not given');
    error('Expected 3D points for plane fitting');
end

% Fill in useful constants
NpointsInFace = length(patchPoints(:,1));
NtestPoints = length(testPoints(:,1));


% Fit a plane to the face
% The plane has the form: 
%             A*(x-x0) + B(y-y0) + C(z-z0) = 0
%
% where A, B, C, are given by the unit_normal_vector: the unit vector in
% the direction of <A, B, C>
% and x0, y0, z0 are given by the base_point

[normalToPlane, orthoBasisForPlane, base_point] = fcn_INTERNAL_affineFit(patchPoints);

% Check ordering of orthoBasisForPlane - not sure if it is XY or YX
crossProduct = cross(orthoBasisForPlane(:,1)',orthoBasisForPlane(:,2)');
dotProduct = sum(crossProduct.*normalToPlane',2);
if dotProduct<0
    orthoBasisForPlane = orthoBasisForPlane(:,[2 1]);
end

magnitude_normal = sum(normalToPlane'.^2,2).^0.5;
unit_normal_vector = normalToPlane'./magnitude_normal;
plane_distances = sum(ones(NpointsInFace,1)*unit_normal_vector.*(patchPoints - ones(NpointsInFace,1)*base_point),2);
standard_deviation_in_plane_orthogonals = std(plane_distances);

if standard_deviation_in_plane_orthogonals>tolerance
    warning('on','backtrace');
    warning('A 3d face given where all points were not in a plane');
    error('Face given that does not have all points in a plane');
end

% Check the height of the testPoints points relative to the plane
testPoint_offPlaneHeight = sum((ones(NtestPoints,1)*unit_normal_vector).*(testPoints - ones(NtestPoints,1)*base_point),2);
flag_isOn3DPatchPlane = abs(testPoint_offPlaneHeight)<tolerance;


% Use the ortho basis vectors to convert testPoints into Cartesian 2D
localXbasis = orthoBasisForPlane(:,1)';
localYbasis = orthoBasisForPlane(:,2)';
localZbasis = unit_normal_vector(1,:);
localBasisVectors = [localXbasis; localYbasis; localZbasis];

localXY_base_point  = [sum(base_point.*localXbasis,2) sum(base_point.*localYbasis,2)];
localXZ_base_point  = [sum(base_point.*localXbasis,2) sum(base_point.*localZbasis,2)];

localXY_patch      = [sum(patchPoints.*localXbasis,2) sum(patchPoints.*localYbasis,2)];
localXZ_patch      = [sum(patchPoints.*localXbasis,2) sum(patchPoints.*localZbasis,2)];
localXY_testPoints = [sum(testPoints.*localXbasis,2)  sum(testPoints.*localYbasis,2)];

localXY_poly = polyshape(localXY_patch);
flag_projectsInsidePatch = isinterior(localXY_poly, localXY_testPoints);

% Set the output
flag_isInside3Dpatch = flag_isOn3DPatchPlane & flag_projectsInsidePatch;

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

    flag_labelPoints = 0; % A flag, for debugging, to label the points

    % Plot the results in point space
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end      
    
 
    % In subplot 1, plot the flane fit in 3D
    subplot(1,3,1);
    hold on;
    grid on;
    % axis equal;
    title('Plane fitting XYZ') 

    % Plot the patchPoints
    plot3(patchPoints(:,1),patchPoints(:,2),patchPoints(:,3),'b.','MarkerSize',10,'DisplayName','patchPoints');

    % Plot the patchPoints as a patch

    patchPointStructure.FaceColor = 'flat';
    patchPointStructure.FaceVertexCData = [0 0 1];
    patchPointStructure.FaceAlpha = 0.3;

    patchPointStructure.Vertices = patchPoints;
    patchPointStructure.Faces = (1:NpointsInFace);
    patch(patchPointStructure,'DisplayName','patchPoints region');

    % Plot the testPoints
    plot3(testPoints(:,1),testPoints(:,2),testPoints(:,3),'k.','MarkerSize',20,'DisplayName','testPoints');

    % Circle planar points and insidePatch points
    plot3(testPoints(flag_isOn3DPatchPlane,1),testPoints(flag_isOn3DPatchPlane,2),testPoints(flag_isOn3DPatchPlane,3), ...
        'm.','MarkerSize',18,'DisplayName','flag_isOn3DPatchPlane');
    plot3(testPoints(flag_projectsInsidePatch,1),testPoints(flag_projectsInsidePatch,2),testPoints(flag_projectsInsidePatch,3), ...
        'c.','MarkerSize',13,'DisplayName','flag_projectsInsidePatch');

    % Plot final results
    pointsToPlot = find(flag_isInside3Dpatch);
    plot3(testPoints(pointsToPlot,1),testPoints(pointsToPlot,2),testPoints(pointsToPlot,3), ...
        'g.','MarkerSize',10,'DisplayName','flag_isInside3Dpatch');

    legend('Interpreter','none');

    view(3);

    % Make axis slightly larger?
    if flag_rescale_axis
        axis equal
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        axis_range_z = temp_axis(6)-temp_axis(5);
        percent_larger = 0.3;
        axis([...
            temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,...
            temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y,...
            temp_axis(5)-percent_larger*axis_range_z, temp_axis(6)+percent_larger*axis_range_z
            ]);
      
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    ax = gca;
    ax.Clipping = "off";

    % Plot the plane
    reshaped_axis = reshape(temp_axis',2,[]);
    maxDifference = sum(diff(reshaped_axis).^2,2).^0.5;
    nudge = maxDifference*0.01;

    % Plot the unit vectors for local basis
    quiver3(ones(3,1)*base_point(1,1),ones(3,1)*base_point(1,2),ones(3,1)*base_point(1,3), ...
        localBasisVectors(:,1),localBasisVectors(:,2),localBasisVectors(:,3),...
        0,'g','Linewidth',2,'DisplayName','local basis vectors');
    text(base_point(1,1)+localBasisVectors(1,1), base_point(1,2)+localBasisVectors(1,2), base_point(1,3)+localBasisVectors(1,3),'x','Color',[0 0.7 0]);
    text(base_point(1,1)+localBasisVectors(2,1), base_point(1,2)+localBasisVectors(2,2), base_point(1,3)+localBasisVectors(2,3),'y','Color',[0 0.7 0]);
    text(base_point(1,1)+localBasisVectors(3,1), base_point(1,2)+localBasisVectors(3,2), base_point(1,3)+localBasisVectors(3,3),'z','Color',[0 0.7 0]);

    % Label the testPoints with their numbers?
    if 1 == flag_labelPoints
        for ith_vertex = 1:length(testPoints(:,1))
            if dimension_of_points==3
                text(testPoints(ith_vertex,1)+nudge,testPoints(ith_vertex,2),testPoints(ith_vertex,3),...
                    sprintf('%.0d',ith_vertex),'Color',[0 0 0]);

            else
                text(testPoints(ith_vertex,1)+nudge,testPoints(ith_vertex,2),...
                    sprintf('%.0d',ith_vertex),'Color',[0 0 0]);
            end

        end
    end


    %% In subplot 2, plot the fit in 2D in local coordinates
    subplot(1,3,2);
    hold on;
    grid on;
    % axis equal;
    title('Local XY plane') 

    % Plot the local_poly
    patchPointStructure.FaceColor = 'flat';
    patchPointStructure.FaceVertexCData = [0 0 1];
    patchPointStructure.FaceAlpha = 0.3;
    patchPointStructure.Vertices = localXY_patch;
    patchPointStructure.Faces = (1:NpointsInFace);
    patch(patchPointStructure,'DisplayName','patchPoints region');


    % % Plot the patchPoints as a patch
    % patchPointStructure.FaceColor = 'flat';
    % patchPointStructure.FaceVertexCData = [0 0 1];
    % patchPointStructure.FaceAlpha = 0.3;
    % 
    % patchPointStructure.Vertices = patchPoints;
    % patchPointStructure.Faces = (1:NpointsInFace);
    % patch(patchPointStructure,'DisplayName','patchPoints region');

    % Plot the testPoints
    plot(localXY_testPoints(:,1),localXY_testPoints(:,2),'k.','MarkerSize',20,'DisplayName','testPoints (local)');

    % Circle planar points and insidePatch points
    plot(localXY_testPoints(flag_isOn3DPatchPlane,1),localXY_testPoints(flag_isOn3DPatchPlane,2), ...
        'm.','MarkerSize',18,'DisplayName','flag_isOn3DPatchPlane (local)');
    plot(localXY_testPoints(flag_projectsInsidePatch,1),localXY_testPoints(flag_projectsInsidePatch,2), ...
        'c.','MarkerSize',13,'DisplayName','flag_projectsInsidePatch (local)');

    % Plot final results
    pointsToPlot = find(flag_isInside3Dpatch);
    plot(localXY_testPoints(pointsToPlot,1),localXY_testPoints(pointsToPlot,2), ...
        'g.','MarkerSize',10,'DisplayName','flag_isInside3Dpatch (local)');

    legend('Interpreter','none');
   
    view(2);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.1;
        axis([...
            temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,...
            temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y...
            ]);
      
    end

    xlabel('local X');
    ylabel('local Y');
    ax = gca;
    ax.Clipping = "off";

    % Plot the plane
    reshaped_axis = reshape(temp_axis',2,[]);
    maxDifference = sum(diff(reshaped_axis).^2,2).^0.5;
    nudge = maxDifference*0.01;

    % Plot the unit vectors for local basis
    quiver(ones(2,1)*localXY_base_point(1,1),ones(2,1)*localXY_base_point(1,2), ...
        [1; 0],[0; 1],...
        0,'g','Linewidth',2,'DisplayName','local basis vectors');
    text(localXY_base_point(1,1)+1, localXY_base_point(1,2)+0,'x','Color',[0 0.7 0]);
    text(localXY_base_point(1,1)+0, localXY_base_point(1,2)+1,'y','Color',[0 0.7 0]);

    % Label the testPoints with their numbers?
    if 1 == flag_labelPoints
        for ith_vertex = 1:length(localXY_testPoints(:,1))
            text(localXY_testPoints(ith_vertex,1)+nudge,localXY_testPoints(ith_vertex,2),...
                sprintf('%.0d',ith_vertex),'Color',[0 0 0]);
        end
    end

    axis equal

    %% In subplot 3, plot the fit in 2D in local coordinates
    subplot(1,3,3);
    hold on;
    grid on;
    % axis equal;
    title('Local XZ plane') 

    % Plot the local_poly
    patchPointStructure.FaceColor = 'flat';
    patchPointStructure.FaceVertexCData = [0 0 1];
    patchPointStructure.FaceAlpha = 0.3;
    patchPointStructure.Vertices = localXZ_patch;
    patchPointStructure.Faces = (1:NpointsInFace);
    patch(patchPointStructure,'DisplayName','patchPoints region');

    % % Plot the patchPoints as a patch
    % patchPointStructure.FaceColor = 'flat';
    % patchPointStructure.FaceVertexCData = [0 0 1];
    % patchPointStructure.FaceAlpha = 0.3;
    % 
    % patchPointStructure.Vertices = patchPoints;
    % patchPointStructure.Faces = (1:NpointsInFace);
    % patch(patchPointStructure,'DisplayName','patchPoints region');

    % Plot the testPoints
    plot(localXY_testPoints(:,1),testPoint_offPlaneHeight,'k.','MarkerSize',20,'DisplayName','testPoints (local)');

    % Circle planar points and insidePatch points
    plot(localXY_testPoints(flag_isOn3DPatchPlane,1),testPoint_offPlaneHeight(flag_isOn3DPatchPlane,1), ...
        'm.','MarkerSize',18,'DisplayName','flag_isOn3DPatchPlane (local)');
    plot(localXY_testPoints(flag_projectsInsidePatch,1),testPoint_offPlaneHeight(flag_projectsInsidePatch,1), ...
        'c.','MarkerSize',13,'DisplayName','flag_projectsInsidePatch (local)');

    % Plot final results
    pointsToPlot = find(flag_isInside3Dpatch);
    plot(localXY_testPoints(pointsToPlot,1),testPoint_offPlaneHeight(pointsToPlot,1), ...
        'g.','MarkerSize',10,'DisplayName','flag_isInside3Dpatch (local)');

    legend('Interpreter','none');
   
    view(2);

    % Make axis slightly larger?
    if flag_rescale_axis
        temp_axis = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp_axis(2)-temp_axis(1);
        axis_range_y = temp_axis(4)-temp_axis(3);
        percent_larger = 0.1;
        axis([...
            temp_axis(1)-percent_larger*axis_range_x, temp_axis(2)+percent_larger*axis_range_x,...
            temp_axis(3)-percent_larger*axis_range_y, temp_axis(4)+percent_larger*axis_range_y...
            ]);
      
    end

    xlabel('local X');
    ylabel('local Z');
    ax = gca;
    ax.Clipping = "off";

    % Plot the plane
    reshaped_axis = reshape(temp_axis',2,[]);
    maxDifference = sum(diff(reshaped_axis).^2,2).^0.5;
    nudge = maxDifference*0.01;

    % Plot the unit vectors for local basis
    quiver(ones(2,1)*localXZ_base_point(1,1),ones(2,1)*localXZ_base_point(1,2), ...
        [1; 0],[0; 1],...
        0,'g','Linewidth',2,'DisplayName','local basis vectors');
    text(localXZ_base_point(1,1)+1, localXZ_base_point(1,2)+0,'x','Color',[0 0.7 0]);
    text(localXZ_base_point(1,1)+0, localXZ_base_point(1,2)+1,'z','Color',[0 0.7 0]);

    % Label the testPoints with their numbers?
    if 1 == flag_labelPoints
        for ith_vertex = 1:length(localXY_testPoints(:,1))
            text(localXY_testPoints(ith_vertex,1)+nudge,testPoint_offPlaneHeight(ith_vertex,1),...
                sprintf('%.0d',ith_vertex),'Color',[0 0 0]);
        end
    end

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