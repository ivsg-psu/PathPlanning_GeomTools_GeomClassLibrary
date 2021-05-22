% script_test_fcn_geometry_findAngleUsing3PointsOnCircle
% Tests fcn_geometry_findAngleUsing3PointsOnCircle

% Revision history:
%      2021_04_28:
%      -- Added more comments, made testing cases more clear


%% BASIC example for one circle, incoming and outgoing are 90 degrees
fig_num = 1;
apex_points = [1 0];
centers = [0 0];
radii = [1]; %#ok<*NBRAK>
start_angles = [45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-45]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [0 2^0.5];
outgoing_destination_points = [0 -2^0.5];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num) %#ok<*ASGLU,*NOPTS>

%% BASIC example for one circle, incoming and outgoing are 270 degrees
fig_num = 2;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 270 degrees
% BUT the apex angle is in the wrong location (so BAD)
fig_num = 3;
apex_points = [-1 0];
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 270 degrees
% BUT the apex angle is in the wierd location
fig_num = 4;
apex_points = [-1 1]/2^0.5;
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 90 degrees
fig_num = 5;
apex_points = [1 1]/(2^0.5);
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [1 -1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)


%% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a good situation
fig_num = 6;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [180]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 -1];
outgoing_destination_points = [1 -1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a good situation
fig_num = 7;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-90]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [-1 -1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 0 degrees,
% and it's a BAD situation
fig_num = 8;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-90]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [1 -1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)


%% BASIC example for one circle, incoming and outgoing are 0 degrees,
% and it's a good situation (grazing contact)
fig_num = 9;
apex_points = [1 0];
centers = [0 0];
radii = [1];
% start_angles = [180]*pi/180;
start_points_on_circle = [1 0];
end_points_on_circle = [1 0];
incoming_source_points = [1 1];
outgoing_destination_points = [1 -1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a BAD situation (grazing contact)
fig_num = 10;
apex_points = [1 0];
centers = [0 0];
radii = [1];
% start_angles = [180]*pi/180;
start_points_on_circle = [1 0];
end_points_on_circle = [1 0];
incoming_source_points = [1 1];
outgoing_destination_points = [1 1];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)

%% BASIC example for one circle that is NOT feasible
fig_num = 11;
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-45]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [0 2^0.5];
outgoing_destination_points = [2^0.5 0];

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    apex_points,...
    centers,...
    start_points_on_circle,...
    end_points_on_circle,...
    radii,...
    incoming_source_points,...
    outgoing_destination_points,fig_num)


%% Fill in vectorized cases 
% Initialize variables
all_apex_points = [];
all_centers = [];
all_radii = [];
all_start_points_on_circle = [];
all_end_points_on_circle = [];
all_incoming_source_points = [];
all_outgoing_destination_points = [];
index = 0;
offset = 4;



apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-45]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [0 2^0.5];
outgoing_destination_points = [0 -2^0.5];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% CASE 1
% BASIC example for one circle, incoming and outgoing are 270 degrees
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% CASE 2
% BASIC example for one circle, incoming and outgoing are 270 degrees
% BUT the apex angle is in the wrong location (so BAD)
apex_points = [-1 0];
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 270 degrees
% BUT the apex angle is in the weird location
apex_points = [-1 1]/2^0.5;
centers = [0 0];
radii = [1];
start_angles = [135]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-135]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-2^0.5 0];
outgoing_destination_points = [-2^0.5 0];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 90 degrees
apex_points = [1 1]/(2^0.5);
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [1 -1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a good situation
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [180]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [0]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 -1];
outgoing_destination_points = [1 -1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a good situation
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-90]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [-1 -1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 0 degrees,
% and it's a BAD situation
apex_points = [1 0];
centers = [0 0];
radii = [1];
start_angles = [90]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-90]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [-1 1];
outgoing_destination_points = [1 -1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 0 degrees,
% and it's a good situation (grazing contact)
apex_points = [1 0];
centers = [0 0];
radii = [1];
% start_angles = [180]*pi/180;
start_points_on_circle = [1 0];
end_points_on_circle = [1 0];
incoming_source_points = [1 1];
outgoing_destination_points = [1 -1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle, incoming and outgoing are 180 degrees,
% and it's a BAD situation (grazing contact)
apex_points = [1 0];
centers = [0 0];
radii = [1];
% start_angles = [180]*pi/180;
start_points_on_circle = [1 0];
end_points_on_circle = [1 0];
incoming_source_points = [1 1];
outgoing_destination_points = [1 1];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

% BASIC example for one circle that is NOT feasible
centers = [0 0];
radii = [1];

apex_angle = 160*pi/180;
apex_points = [radii.*cos(apex_angle) radii.*sin(apex_angle)];
start_angles = [45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-45]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
incoming_source_points = [0 2^0.5];
outgoing_destination_points = [2^0.5 0];

shift = [0 index*offset];
all_apex_points = [all_apex_points; apex_points+shift];
all_centers = [all_centers; centers+shift];
all_radii = [all_radii; radii];
all_start_points_on_circle = [all_start_points_on_circle; start_points_on_circle+shift];
all_end_points_on_circle = [all_end_points_on_circle; end_points_on_circle+shift];
all_incoming_source_points = [all_incoming_source_points; incoming_source_points+shift];
all_outgoing_destination_points = [all_outgoing_destination_points; outgoing_destination_points+shift];
index = index+1;

[angles, better_angles] = fcn_geometry_findAngleUsing3PointsOnCircle(...
    all_apex_points,...
    all_centers,...
    all_start_points_on_circle,...
    all_end_points_on_circle,...
    all_radii,...
    all_incoming_source_points,...
    all_outgoing_destination_points,23);


