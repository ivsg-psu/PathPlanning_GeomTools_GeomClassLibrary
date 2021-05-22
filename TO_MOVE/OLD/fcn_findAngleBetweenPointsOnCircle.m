function [angles, better_angles,better_angle_range,inpoints_are_closer_to_apex] = fcn_findAngleBetweenPointsOnCircle(apex_points, centers,start_points_on_circle,end_points_on_circle,radii,incoming_source_points,outgoing_destination_points,varargin)
% fcn_findAngleBetweenPointsOnCircle - finds the angle between two points
% on a circle
% angles = fcn_findAngleBetweenPointsOnCircle(apex_points, centers,start_points_on_circle,end_points_on_circle,radii,incoming_source_points,outgoing_destination_points,varargin)
% This function calculates the angle from the start_points location to the
% end_points, in the direction of the vector given by is_clockwise. The
% inputs are:
% apex_points = an Nx2 vector in [x y] format of the apex points
% centers = an Nx2 vector in [x y] of the points of circle centers
% start_points_on_circle = an Nx2 vector in [x y] of the points where sectors start
% end_points_on_circle = an Nx2 vector in [x y] of the points where sectors end
% radii = an Nx1 vector of the radii of the circles (to avoid calculation
% time)
% incoming_source_points = an Nx2 vector in [x y] of the points where the
% incoming line is originating from
% outgoing_destination_points = an Nx2 vector in [x y] of the points where
% outgoing line segment is going to.
%
% Examples:
%
%      % BASIC example for one circle, incoming and outgoing are 90 degrees
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [45]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-45]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [0 2^0.5];
%      outgoing_destination_points = [0 -2^0.5];
%
%      % BASIC example for one circle, incoming and outgoing are 270 degrees
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [135]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-135]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-2^0.5 0];
%      outgoing_destination_points = [-2^0.5 0];
%
%      % BASIC example for one circle, incoming and outgoing are 270 degrees
%      % BUT the apex angle is in the wrong location (so BAD)
%      apex_points = [-1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [135]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-135]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-2^0.5 0];
%      outgoing_destination_points = [-2^0.5 0];
%
%      % BASIC example for one circle, incoming and outgoing are 270 degrees
%      % BUT the apex angle is in the wierd location 
%      apex_points = [-1 1]/2^0.5;
%      centers = [0 0];
%      radii = [1];
%      start_angles = [135]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-135]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-2^0.5 0];
%      outgoing_destination_points = [-2^0.5 0];
%
%      % BASIC example for one circle, incoming and outgoing are 90 degrees
%      apex_points = [1 1]/(2^0.5);
%      centers = [0 0];
%      radii = [1];
%      start_angles = [90]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [0]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-1 1];
%      outgoing_destination_points = [1 -1];
%
%      % BASIC example for one circle, incoming and outgoing are 180 degrees,
%      % and it's a good situation
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [180]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [0]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-1 -1];
%      outgoing_destination_points = [1 -1];
%
%      % BASIC example for one circle, incoming and outgoing are 180 degrees,
%      % and it's a good situation
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [90]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-90]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-1 1];
%      outgoing_destination_points = [-1 -1];
%
%      % BASIC example for one circle, incoming and outgoing are 0 degrees,
%      % and it's a BAD situation
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [90]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-90]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [-1 1];
%      outgoing_destination_points = [1 -1];
%
%
%      % BASIC example for one circle, incoming and outgoing are 0 degrees,
%      % and it's a good situation (grazing contact)
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [180]*pi/180;
%      start_points_on_circle = [1 0];
%      end_points_on_circle = [1 0];
%      incoming_source_points = [1 1];
%      outgoing_destination_points = [1 -1];
%
%      % BASIC example for one circle, incoming and outgoing are 180 degrees,
%      % and it's a BAD situation (grazing contact)
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [180]*pi/180;
%      start_points_on_circle = [1 0];
%      end_points_on_circle = [1 0];
%      incoming_source_points = [1 1];
%      outgoing_destination_points = [1 1];
%
%      % BASIC example for one circle that is NOT feasible
%      apex_points = [1 0];
%      centers = [0 0];
%      radii = [1];
%      start_angles = [45]*pi/180;
%      start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
%      end_angles = [-45]*pi/180;
%      end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
%      incoming_source_points = [0 2^0.5];
%      outgoing_destination_points = [2^0.5 0];
%
%      angles = fcn_findAngleBetweenPointsOnCircle(centers,start_points_on_circle,end_points_on_circle,radii,incoming_source_points,outgoing_destination_points,varargin)
%
% This function was written on 2020_03_28 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%
% See: http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm for
% the mathematics being implemented here.

do_debug = 0;

%% check input arguments to see if we need to make a figure
if nargin < 7 || nargin > 8
    error('Incorrect number of input arguments')
end

if 8 == nargin % Did user give a figure number?
    fig_num = varargin{1};
    figure(fig_num);
    do_debug = 1; % To force plotting
else
    if do_debug  % Only create a new figure if debugging
        fig = figure(1); %#ok<UNRCH> % create new figure with next default index
        fig_num = get(fig,'Number');
    end
end

%% Set up variables
num_circles = length(radii(:,1));

%% Error check all the inputs
% Check the arguments themselves
% Check if all are in correct column format
if length(apex_points(1,:))~=2
    error('The apex_points vector is expected to have 2 columns, specifying x and y locations of the apexes respectively');
end
if length(centers(1,:))~=2
    error('The centers is expected to have 2 columns, specifying x and y locations of the pivot circle centers respectively');
end
if length(start_points_on_circle(1,:))~=2
    error('The starting circles is expected to have 2 columns, specifying x and y locations of the centers respectively');
end
if length(end_points_on_circle(1,:))~=2
    error('The ending circles is expected to have 2 columns, specifying x and y locations of the centers respectively');
end
if length(radii(1,:))~=1
    error('The radii vector for the pivot circles should have only 1 column, specifying radii for each circle respectively');
end
if length(incoming_source_points(1,:))~=2
    error('The incoming_source_points vector is expected to have 2 columns, specifying x and y start locations of the incoming line segment.');
end
if length(outgoing_destination_points(1,:))~=2
    error('The outgoing_destination_points vector is expected to have 2 columns, specifying x and y final locations of the outgoing line segment.');
end

% Check the lengths of each
if length(apex_points(:,1))~= num_circles
    error('The number of apex points must match the number of pivot circles.');
end
if length(centers(:,1))~= num_circles
    error('The number of radii must match the number of pivot circles.');
end
if (length(start_points_on_circle(:,1))~=num_circles) || (length(end_points_on_circle(:,1))~=num_circles)
    error('The column length of the number of start/end points on the circle must match the number of circles.');
end
if (length(incoming_source_points(:,1))~=num_circles) || (length(outgoing_destination_points(:,1))~=num_circles)
    error('The column length of the number of incoming and outgoing points for each circle must match the number of circles.');
end

% If made it here, the vectors are probably good!

%% Step 0: calculate unit vectors for incoming and outgoing vectors
% [angles, better_angles,better_angle_range,inpoints_are_closer_to_apex] = 
% fcn_findAngleBetweenPointsOnCircle(apex_points, centers,start_points_on_circle,...
% end_points_on_circle,radii,incoming_source_points,outgoing_destination_points,varargin)
% Makes things easy later
diff_in = start_points_on_circle - incoming_source_points;
diff_out = outgoing_destination_points - end_points_on_circle;
mag_in = sum(diff_in.^2,2).^0.5;
mag_out = sum(diff_out.^2,2).^0.5;
unit_vector_in = diff_in./mag_in;
unit_vector_out = diff_out./mag_out;
unit_radial_to_inpoints = (start_points_on_circle - centers)./radii;
unit_radial_to_outpoints = (end_points_on_circle - centers)./radii;
mean_unit_vector_in_and_out = (unit_radial_to_inpoints+unit_radial_to_outpoints)./2;
mean_unit_vector_in_and_out = mean_unit_vector_in_and_out./sum(mean_unit_vector_in_and_out.^2,2);
    
%% Plot situation at start? (for deep debugging!)
if 1==0
    figure(1111);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    plot(apex_points(:,1),apex_points(:,2),'r*');
    plot_circle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    plot(start_points_on_circle(:,1),start_points_on_circle(:,2),'go');
    plot(end_points_on_circle(:,1),end_points_on_circle(:,2),'ro');
    
    
    for i=1:num_circles
        plot([incoming_source_points(i,1); start_points_on_circle(i,1)],...
            [incoming_source_points(i,2); start_points_on_circle(i,2)],...
            'g-');
        plot([end_points_on_circle(i,1); outgoing_destination_points(i,1)],...
            [end_points_on_circle(i,2);  outgoing_destination_points(i,2)],...
            'r-');
        
        % Plot unit vectors
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_inpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_inpoints(i,2)],'c');
        
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_outpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_outpoints(i,2)],'m');
        
        plot([centers(i,1); centers(i,1)+radii(i).*mean_unit_vector_in_and_out(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*mean_unit_vector_in_and_out(i,2)],'g');

    end
    
    
end



%% Step 1: calculate the dot product angle from in and out unit vectors
dot_product = sum(unit_vector_in.*unit_vector_out,2);
angles = acos(dot_product);

%% Step 2: check weird situations
% Check 1: are they parallel and in same direction? (not possible)
% If in the same direction, that's not feasible unless the apex point, in
% point, and outpout are all the same point (e.g. a grazing contact)
vectors_are_same = unit_vector_in==unit_vector_out;
vectors_are_same = vectors_are_same(:,1).*vectors_are_same(:,2);
indices = find(vectors_are_same==1);
for i=1:length(indices)
    current_index = indices(i);
    if ~isequal(apex_points(current_index,:),start_points_on_circle(current_index,:)) && ~isequal(apex_points(current_index,:),end_points_on_circle(current_index,:))
        angles(current_index) = NaN; % This is a bad situation
    end
end

% Check 2: they are 180 degrees apart, and go in/out of same point
indices_180_apart = find(dot_product==-1);
for i=1:length(indices_180_apart)
    current_index = indices_180_apart(i);
    if isequal(start_points_on_circle(current_index,:),end_points_on_circle(current_index,:))
        angles(current_index) = NaN; % This is a bad situation
    end
end

% Check 3: Check to see if vectors "spin" the circle in the same direction
% Do the cross product of the radial vectors with the incoming and outgoing
% vectors

cross1 = cross([unit_radial_to_inpoints, zeros(num_circles,1)],...
    [unit_vector_in, zeros(num_circles,1)]);
cross_radius_out_to_out = cross([unit_radial_to_outpoints, zeros(num_circles,1)],...
    [unit_vector_out, zeros(num_circles,1)]);

% If the cross products have different signs, then they are bad
bad_cross_product_indices = (cross1(:,3).*cross_radius_out_to_out(:,3))< 0;
angles(bad_cross_product_indices) = NaN;

%% Step 3: fix angles for situations where the sector is larger than 180
% The angle we are looking for is the reflex angle, which we can find by
% checking the directions of the cross products
cross_in_to_radius_in = cross(...
    [unit_vector_in, zeros(num_circles,1)],...
    [unit_radial_to_inpoints, zeros(num_circles,1)]);
cross_in_to_out = cross(...
    [unit_vector_in, zeros(num_circles,1)],...
    [unit_vector_out, zeros(num_circles,1)]);


% If the cross_in_to_radius_in is in opposite direction from
% cross_radius_to_out, then the angles we are calculating are incorrect. We
% need to take the refelex angle.
need_reflex_angles = (cross_in_to_radius_in(:,3)...
    .*cross_in_to_out(:,3))>0;
angles(need_reflex_angles)=2*pi - angles(need_reflex_angles);


%% Step 4: confirm that the apex angle is in a valid part of the sector
% Calculate unit vectors from entry point in circle to exit point in circle
% and unit vector from entry point to the apex
diff_in_to_out  = end_points_on_circle - start_points_on_circle;
diff_in_to_apex = apex_points - start_points_on_circle;
mag_in_to_out = sum(diff_in_to_out.^2,2).^0.5;
mag_in_to_apex = sum(diff_in_to_apex.^2,2).^0.5;
unit_vector_in_to_out = diff_in_to_out./mag_in_to_out;
unit_vector_in_to_apex = diff_in_to_apex./mag_in_to_apex;

dot_product_in_to_in_to_out = sum(unit_vector_in.*unit_vector_in_to_out,2);
dot_product_in_to_in_to_apex = sum(unit_vector_in.*unit_vector_in_to_apex,2);
angles1 = acos(dot_product_in_to_in_to_out);
angles2 = acos(dot_product_in_to_in_to_apex);

% If angle to apex point larger than angle to exit point, then it's not in
% the arc
angles(angles2>angles1)=NaN;

% If the angles is bad, then give a better option. First, fill in
% better_angle suggestions with what the user gave, e.g. apex angles:
vector_center_to_apex = apex_points - centers;
unit_radial_to_apex = vector_center_to_apex./hypot(vector_center_to_apex(:,1),vector_center_to_apex(:,2));
current_angle_apex = mod(pi+atan2(unit_radial_to_apex(:,2),unit_radial_to_apex(:,1)),2*pi);

better_angles = current_angle_apex;
angle_in  = mod(pi+atan2(unit_radial_to_inpoints(:,2),unit_radial_to_inpoints(:,1)),2*pi);
angle_out = mod(pi+atan2(unit_radial_to_inpoints(:,2),unit_radial_to_outpoints(:,1)),2*pi);
better_angle_range = sort([angle_in angle_out],2);
inpoints_are_closer_to_apex = ones(length(centers(:,1)),1);

if any(isnan(angles))  % calculate a better suggestion
    % OLD METHOD: find closest angle
    %     distance_apex_to_inpoints = sum((unit_radial_to_apex - unit_radial_to_inpoints).^2,2);
    %     distance_apex_to_outpoints = sum((unit_radial_to_apex - unit_radial_to_outpoints).^2,2);
    %     in_is_closest = distance_apex_to_inpoints<distance_apex_to_outpoints;
    %     angle_in = atan2(unit_radial_to_inpoints(:,2),unit_radial_to_inpoints(:,1));
    %     angle_out = atan2(unit_radial_to_inpoints(:,2),unit_radial_to_outpoints(:,1));
    %     better_angles(angles2>angles1) = mod(pi+(angle_in(angles2>angles1).*in_is_closest(angles2>angles1) + angle_out(angles2>angles1).*(1-in_is_closest(angles2>angles1))),2*pi);

    % NEW METHOD: take mean of unit vectors to create a new unit vector        
    mean_angle = atan2(mean_unit_vector_in_and_out(:,2),mean_unit_vector_in_and_out(:,1));
    better_angles(angles2>angles1) = mod(pi+(mean_angle(angles2>angles1,:)),2*pi);
    better_angles(bad_cross_product_indices) = NaN; % Do not suggest better if cross product is bad!

    % Calculate distance between apex to inpoints and outpoints
    distance_apex_to_inpoints = sum((unit_radial_to_apex - unit_radial_to_inpoints).^2,2);
    distance_apex_to_outpoints = sum((unit_radial_to_apex - unit_radial_to_outpoints).^2,2);
    inpoints_are_closer_to_apex = 2*(distance_apex_to_inpoints<distance_apex_to_outpoints) - 1;

end

% NOTE: the better angles calculation is giving the wrong sign. Need to
% trace where this is coming from!

%% Plot results?
if do_debug
    figure(fig_num);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    plot(apex_points(:,1),apex_points(:,2),'r*');
    plot_circle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    plot(start_points_on_circle(:,1),start_points_on_circle(:,2),'go');
    plot(end_points_on_circle(:,1),end_points_on_circle(:,2),'go');
    
    for i=1:num_circles
        plot([incoming_source_points(i,1); start_points_on_circle(i,1)],...
            [incoming_source_points(i,2); start_points_on_circle(i,2)],...
            'g-');
        plot([end_points_on_circle(i,1); outgoing_destination_points(i,1)],...
            [end_points_on_circle(i,2);  outgoing_destination_points(i,2)],...
            'r-');

        % Plot unit vectors
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_inpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_inpoints(i,2)],'c');
        
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_outpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_outpoints(i,2)],'m');
        
        plot([centers(i,1); centers(i,1)+radii(i).*mean_unit_vector_in_and_out(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*mean_unit_vector_in_and_out(i,2)],'g');
        
        location = apex_points(i,:);
        text(location(1,1),location(1,2),sprintf(':  %.1f deg',angles(i,1)*180/pi));
    end
end


end

function plot_circle(centers,radii)    % plot all the circle fits
angles = 0:0.01:2*pi;
for i_fit = 1:length(centers(:,1))
    x_circle = centers(i_fit,1) + radii(i_fit) * cos(angles);
    y_circle = centers(i_fit,2) + radii(i_fit) * sin(angles);
    plot(x_circle,y_circle,'b-');
end
end

