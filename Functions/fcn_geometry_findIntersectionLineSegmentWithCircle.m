function [intAngle,intPoint] = fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R)
% fcn_geometry_findIntersectionLineSegmentWithCircle Evaluates an edge
% against a circle to determine whether there will be any intersections
% between the two.
%
% FORMAT:
%
%       [intAngle,intPoint] = fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R)
%
% INPUTS:
%
%      pa: a 1 x 2 vector containing the (x,y) coordinates of one end of
%           the line segment
%      pb: a 1 x 2 vector containing the (x,y) coordinates of one end of
%           the line segment
%      pc: a 1 x 2 vector containing the (x,y) coordinates of the center of
%           the circle
%      R: a scalar parameter containing the radius of the circle
%
% OUTPUTS:
%
%      intAngle: an N x 1 vector of intersection angles, relative to the
%           positive x-axis, with length 0, 1, or 2, depending on the
%           geometry
%      intPoint: an N x 2 vector of intersection points, with length 0, 1,
%           or 2, depending on the geometry
%
% EXAMPLES:
%
%       See the script: script_test_fcn_geometry_findIntersectionLineSegmentWithCircle.m
%       for a full test suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_03_06
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    
end

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs == 1
    % Are there the right number of inputs?
    narginchk(4,4);

    % Check to see that there were three points provided in the proper
    % format
    if any(size(pa) ~= [1 2]) || any(size(pb) ~= [1 2]) || any(size(pc) ~= [1 2])
        warning('One or more input points is the wrong size (should be 1 x 2).');
        return;
    end
    
    % Check the radius
    if any(size(R) ~= [1 1]) || R <= 0
        error('Radius either negative or not scalar valued.')
    end
end

% Compute the point associated with the min radius. This is computed from
% the derivative of the distance expression with respect to a free
% parameter describing the distance along the line segment from point a to
% point b
alpha = -((pa(1)-pc(1))*(pb(1)-pa(1)) + (pa(2)-pc(2))*(pb(2)-pa(2)))/((pb(1)-pa(1))^2 + (pb(2)-pa(2))^2);
% If the minimum radius lies on the segment between points a and b,
% calculate the location of the point of minimum radius
if alpha < 1 && alpha > 0
    pminR = pa + alpha*(pb - pa);
% If the minimum radius is outside of the segment, set the point to
% infinity so that it is ignored when determining intersections
else
    pminR = [inf inf];
end

% Check for the case where there is a tangency between the line segment and
% the circle
if norm(pminR - pc) == R
        % Calculate the associated circular angle of intersection from the
    % scalar distance along the line segment (alpha)
    intAngle(1) = atan2(pminR(2) - pc(2),pminR(1) - pc(1));
    % Normalize the range of the intersection angle to [0,2*pi]
    while intAngle(1) > 2*pi
        intAngle(1) = intAngle(1) - 2*pi;
    end
    while intAngle(1) < 0
        intAngle(1) = intAngle(1) + 2*pi;
    end
    % Compute the actual point of intersection
    intPoint(1,:) = pminR;
    
elseif (norm(pminR-pc) < R && norm(pa-pc) > R && norm(pb-pc) > R)
    % In this case, there are two intersections. One lies between pa and
    % pminR and the other lies between pminR and pb. Compute both.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First line segment, from pa to pminR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the coefficients of the quadratic equation for alpha
    % (the scalar distance along the line segment)
    quadCoefs(1) = (pminR(1) - pa(1))^2 + (pminR(2) - pa(2))^2;
    quadCoefs(2) = 2*((pa(1)-pc(1))*(pminR(1) - pa(1)) + (pa(2)-pc(2))*(pminR(2) - pa(2)));
    quadCoefs(3) = (pa(1)-pc(1))^2 + (pa(2)-pc(2))^2 - R^2;
    % Calculate the additive root for alpha
    alpha = (-quadCoefs(2) + sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    % Check to see if the solution should come from the subtractive root
    % and adjust if necessary
    if alpha > 1
        alpha = (-quadCoefs(2) - sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    end
    % Confirm that we got a good solution. (This shouldn't fail with
    % the previous logic, but it's here as extra validation.)
    if alpha > 1 || alpha < 0
        error('Calculation error, out of range');
    end
    % Calculate the associated circular angle of intersection from the
    % scalar distance along the line segment (alpha)
    intAngle(1) = atan2(pa(2) - pc(2) + alpha*(pminR(2)-pa(2)),...
        pa(1) - pc(1) + alpha*(pminR(1)-pa(1)));
    % Normalize the range of the intersection angle to [0,2*pi]
    while intAngle(1) > 2*pi
        intAngle(1) = intAngle(1) - 2*pi;
    end
    while intAngle(1) < 0
        intAngle(1) = intAngle(1) + 2*pi;
    end
    % Compute the actual point of intersection
    intPoint(1,:) = pa + alpha*(pminR-pa);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Second line segment, from pminR to pb
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the coefficients of the quadratic equation for alpha
    % (the scalar distance along the line segment)
    quadCoefs(1) = (pb(1) - pminR(1))^2 + (pb(2) - pminR(2))^2;
    quadCoefs(2) = 2*((pminR(1)-pc(1))*(pb(1) - pminR(1)) + (pminR(2)-pc(2))*(pb(2) - pminR(2)));
    quadCoefs(3) = (pminR(1)-pc(1))^2 + (pminR(2)-pc(2))^2 - R^2;
    % Calculate the additive root for alpha
    alpha = (-quadCoefs(2) + sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    % Check to see if the solution should come from the subtractive root
    % and adjust if necessary
    if alpha > 1
        alpha = (-quadCoefs(2) - sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    end
    % Confirm that we got a good solution. (This shouldn't fail with
    % the previous logic, but it's here as extra validation.)
    if alpha > 1 || alpha < 0
        error('Calculation error, out of range');
    end
    % Calculate the associated circular angle of intersection from the
    % scalar distance along the line segment (alpha)
    intAngle(2) = atan2(pminR(2) - pc(2) + alpha*(pb(2)-pminR(2)),...
        pminR(1) - pc(1) + alpha*(pb(1)-pminR(1)));
    % Normalize the range of the intersection angle to [0,2*pi]
    while intAngle(2) > 2*pi
        intAngle(2) = intAngle(2) - 2*pi;
    end
    while intAngle(2) < 0
        intAngle(2) = intAngle(2) + 2*pi;
    end
    % Compute the actual point of intersection
    intPoint(2,:) = pminR + alpha*(pb-pminR);
% Next, check to see if no intersections are possible
elseif (norm(pa-pc) > R && norm(pb-pc) > R) || (norm(pa-pc) < R && norm(pb-pc) < R)
    % If no intersection (or two intersections) possible, return empty
    % vectors of the intersection angle and point
    intAngle = [];
    intPoint = [];
% Next, check for the case where there are two intersections between the
% line segment and circle. This occurs if both of the endpoints of the line
% segment are outside of the circle but the distance to the nearest point
% on the line is less than the circle radius.
else
    % Calculate the coefficients of the quadratic equation for alpha
    % (the scalar distance along the line segment)
    quadCoefs(1) = (pb(1) - pa(1))^2 + (pb(2) - pa(2))^2;
    quadCoefs(2) = 2*((pa(1)-pc(1))*(pb(1) - pa(1)) + (pa(2)-pc(2))*(pb(2) - pa(2)));
    quadCoefs(3) = (pa(1)-pc(1))^2 + (pa(2)-pc(2))^2 - R^2;
    % Calculate the additive root for alpha
    alpha = (-quadCoefs(2) + sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    % Check to see if the solution should come from the subtractive root
    % and adjust if necessary
    if alpha > 1
        alpha = (-quadCoefs(2) - sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    end
    % Confirm that we got a good solution. (This shouldn't fail with
    % the previous logic, but it's here as extra validation.)
    if alpha > 1 || alpha < 0
        error('Calculation error, out of range');
    end
    % Calculate the associated circular angle of intersection from the
    % scalar distance along the line segment (alpha)
    intAngle = atan2(pa(2) - pc(2) + alpha*(pb(2)-pa(2)),...
        pa(1) - pc(1) + alpha*(pb(1)-pa(1)));
    % Normalize the range of the intersection angle to [0,2*pi]
    while intAngle > 2*pi
        intAngle = intAngle - 2*pi;
    end
    while intAngle < 0
        intAngle = intAngle + 2*pi;
    end
    % Compute the actual point of intersection
    intPoint = pa + alpha*(pb-pa);
end
end