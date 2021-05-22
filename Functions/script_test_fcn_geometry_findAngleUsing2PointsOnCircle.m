% script_test_fcn_geometry_findAngleUsing2PointsOnCircle
% Tests fcn_geometry_findAngleUsing2PointsOnCircle

% Revision history:
%      2021_05_22:
%      -- Edited for new function name


%% BASIC example for one circle, incoming and outgoing are 90 degrees
fig_num = 1;
centers = [0 0];
radii = [1]; 

start_angles = [45]*pi/180;
start_points_on_circle = [radii.*cos(start_angles) radii.*sin(start_angles)];
end_angles = [-45]*pi/180;
end_points_on_circle = [radii.*cos(end_angles) radii.*sin(end_angles)];
cross_products = [1];

true_angle = start_angles - end_angles;

[angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products,...
    fig_num);

fcn_summarize(angles,... 
    true_angle,...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products);


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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
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

[angles, better_angles] = fcn_geometry_findAngleUsing2PointsOnCircle(...
    all_apex_points,...
    all_centers,...
    all_start_points_on_circle,...
    all_end_points_on_circle,...
    all_radii,...
    all_incoming_source_points,...
    all_outgoing_destination_points,23);

%% Function start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_summarize(...
    angles,...    
    true_angle,...
    centers,...
    radii,...
    start_points_on_circle,...
    end_points_on_circle,...
    cross_products)

for i=1:length(angles)
    fprintf(1,'True angle: %.2f\n',true_angle(i,1)*180/pi);
    fprintf(1,'Centers: %.2f %.2f\n',centers(i,1),centers(i,1));
    fprintf(1,'Raddi: %.2f ',radii(i,1));
    fprintf(1,'Start points: %.2f %.2f\n',start_points_on_circle(i,1),start_points_on_circle(i,1));
    fprintf(1,'End points: %.2f %.2f\n',end_points_on_circle(i,1),end_points_on_circle(i,1));
    fprintf(1,'Cross products: %.2f \n',cross_products(i,1));
    fprintf(1,'Calculated angle: %.2f\n',angles(i,1)*180/pi);

end % Ends for loop
end % Ends function

