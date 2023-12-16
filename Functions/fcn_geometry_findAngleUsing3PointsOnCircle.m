function [...
    angles,...
    better_angles,...
    better_angle_range,...
    inpoints_are_closer_to_apex] ...
    = ...
    fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,...
    varargin)

% fcn_geometry_findAngleUsing3PointsOnCircle -  This function calculates
% the angle from the start_points location to the end_points, in the
% direction of the vector given by is_clockwise.
%
% FORMAT:
%
% [...
%     angles,...
%     better_angles,...
%     better_angle_range,...
%     inpoints_are_closer_to_apex] ...
%     = ...
%     fcn_geometry_findAngleUsing3PointsOnCircle(...
%     apex_points,...
%     centers,...
%     start_points_on_circle,...
%     end_points_on_circle,...
%     radii,...
%     incoming_source_points,...
%     outgoing_destination_points,...
%     varargin)
%
% INPUTS:
%
%      apex_points: an [N x 2] vector of X,Y data for each apex point
%
%      centers: an [N x 2] vector in [x y] of the points of circle centers
%
%      start_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors start
%
%      end_points_on_circle: an [N x 2] vector in [x y] of the points
%      where sectors end
%
%      radii: a [N x 1] vector of the radii of the circles (to avoid
%      calculation time)
%
%      incoming_source_points: an [N x 2] vector in [x y] of the points
%      where the incoming line is originating from
%
%      outgoing_destination_points: an [N x 2] vector in [x y] of the points
%      where the outgoing line segment is going to.
%
%
%      (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%
% DEPENDENCIES:
%
%      fcn_geometry_checkInputsToFunctions
%      fcn_geometry_findTangentPointFromPointToCircle
%      fcn_geometry_findTangentPointsFromPointToCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_findAngleUsing3PointsOnCircle
% for a full test suite.
%
% This function was written on 2020_03_28 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%
% See: http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm for
% the mathematics being implemented here.

% Revision History:
% 2021-04-22
% -- Added input checking
% -- Added cross-product to see if innner/outer connected


%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

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

if flag_check_inputs    
    % Are there the right number of inputs?
    narginchk(7,8);
    
    % Check the apex_points input
    fcn_geometry_checkInputsToFunctions(...
        apex_points, '2column_of_numbers');
    
    num_circles = length(apex_points(:,1));
    
    % Check the centers input
    fcn_geometry_checkInputsToFunctions(...
        centers, '2column_of_numbers',num_circles);
    
    % Check the start_points_on_circle input
    fcn_geometry_checkInputsToFunctions(...
        start_points_on_circle, '2column_of_numbers',num_circles);
    
    % Check the end_points_on_circle input
    fcn_geometry_checkInputsToFunctions(...
        end_points_on_circle, '2column_of_numbers',num_circles);
    
    % Check the internal_apex_angles input
    fcn_geometry_checkInputsToFunctions(...
        radii, 'column_of_numbers',num_circles);
    
    % Check the incoming_source_points input
    fcn_geometry_checkInputsToFunctions(...
        incoming_source_points, '2column_of_numbers',num_circles);
    
    % Check the outgoing_destination_points input
    fcn_geometry_checkInputsToFunctions(...
        outgoing_destination_points, '2column_of_numbers',num_circles);
    
end
    

% Does user want to show the plots?
if 8 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plot = 1;
    flag_new_figure = 0;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
        flag_new_figure = 1;
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



%% Step 0: calculate unit vectors for incoming and outgoing vectors
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
    
% %% Plot situation at start? (for deep debugging!)
% if 1==0
%     figure(1111);
%     clf;
%     hold on;
%     axis equal;
%     grid on; grid minor;
%     
%     % Plot the apex points
%     plot(apex_points(:,1),apex_points(:,2),'r*');
%     
%     % Plot the circles
%     plot_circle(centers,radii);
%     plot(centers(:,1),centers(:,2),'kx');
%     
%     % Plot the start and end points on the circl
%     plot(start_points_on_circle(:,1),start_points_on_circle(:,2),'go');
%     plot(end_points_on_circle(:,1),end_points_on_circle(:,2),'ro');
%     
%     
%     for i=1:num_circles
%         plot([incoming_source_points(i,1); start_points_on_circle(i,1)],...
%             [incoming_source_points(i,2); start_points_on_circle(i,2)],...
%             'g-');
%         plot([end_points_on_circle(i,1); outgoing_destination_points(i,1)],...
%             [end_points_on_circle(i,2);  outgoing_destination_points(i,2)],...
%             'r-');
%         
%         % Plot unit vectors
%         plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_inpoints(i,1)],...
%             [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_inpoints(i,2)],'c');
%         
%         plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_outpoints(i,1)],...
%             [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_outpoints(i,2)],'m');
%         
%         plot([centers(i,1); centers(i,1)+radii(i).*mean_unit_vector_in_and_out(i,1)],...
%             [centers(i,2); centers(i,2)+radii(i).*mean_unit_vector_in_and_out(i,2)],'g');
% 
%     end
%     
%     
% end



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
    % Set up the figure
    figure(fig_num);
    clf;
    hold on;
    axis equal;
    grid on; grid minor;
    
    
    % Plot the apex points
    plot(apex_points(:,1),apex_points(:,2),'r*');
    
    % Plot the circles
    plot_circle(centers,radii);
    plot(centers(:,1),centers(:,2),'kx');
    
    % Plot the start and end points
    plot(start_points_on_circle(:,1),start_points_on_circle(:,2),'go');
    text(start_points_on_circle(:,1),start_points_on_circle(:,2),'Start');
    
    plot(end_points_on_circle(:,1),end_points_on_circle(:,2),'rx');
    text(end_points_on_circle(:,1),end_points_on_circle(:,2),'End');
    
    for i=1:num_circles
        plot([incoming_source_points(i,1); start_points_on_circle(i,1)],...
            [incoming_source_points(i,2); start_points_on_circle(i,2)],...
            'g-');
        plot([end_points_on_circle(i,1); outgoing_destination_points(i,1)],...
            [end_points_on_circle(i,2);  outgoing_destination_points(i,2)],...
            'r-');

        % Plot unit vectors
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_inpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_inpoints(i,2)],'g');
        
        plot([centers(i,1); centers(i,1)+radii(i).*unit_radial_to_outpoints(i,1)],...
            [centers(i,2); centers(i,2)+radii(i).*unit_radial_to_outpoints(i,2)],'r');
        
        % Plot the angle value
        location = apex_points(i,:);
        text(location(1,1),location(1,2),sprintf(':  %.1f deg',angles(i,1)*180/pi));
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end


end % Ends the function

function plot_circle(centers,radii)    % plot all the circle fits
angles = 0:0.01:2*pi;
for i_fit = 1:length(centers(:,1))
    x_circle = centers(i_fit,1) + radii(i_fit) * cos(angles);
    y_circle = centers(i_fit,2) + radii(i_fit) * sin(angles);
    plot(x_circle,y_circle,'b-');
end
end

