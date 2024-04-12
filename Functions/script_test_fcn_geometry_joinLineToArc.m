%% script_test_fcn_geometry_joinLineToArc
% Exercises the function: fcn_geometry_joinLineToArc
% Revision history:
% 2024_04_12
% -- wrote the code

close all;

%% BASIC test

% Bit 1: line oriented, 
% Bit 2: arc oriented, 
% Bit 3: ends aligned, 
% Bit 4: arc is positive

for ith_test = 1:(2*2*2*2+1)

    numDigits = 4;
    binary_string = dec2bin(ith_test-1,numDigits);

    fig_num = ith_test;
    figure(fig_num); clf;
    title_string = '';

    % URHERE
    if strcmp(binary_string(1),'1')
        line_unit_tangent_vector = [1 0];
        line_base_point_xy       = [0 0];
        line_s_start             = -1;
        line_s_end               = 0;
    end

    if 0==mod(ith_test,2)
        % Arc is positive
        arc_center_xy            = [0 1+error];
        arc_radius               = 1;
        arc_angles               = [270 360]*pi/180;
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 1;
    else
        % Arc is negative
        arc_center_xy            = [0 -1+error];
        arc_radius               = 1;
        arc_angles               = [90 0]*pi/180;
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 0;
    end

    true_line_unit_tangent_vector = [1 0];
    true_start_point_xy = [-1 0];
    true_arc_center_xy  = [0 1];
    true_arc_is_counter_clockwise = 1;

    % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
    line_parameters(1,1:2) = line_unit_tangent_vector;
    line_parameters(1,3:4) = line_base_point_xy;
    line_parameters(1,5)   = line_s_start;
    line_parameters(1,6)   = line_s_end;

    % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
    arc_parameters(1,1:2) = arc_center_xy;
    arc_parameters(1,3)   = arc_radius;
    arc_parameters(1,4:5) = arc_angles;
    arc_parameters(1,7)   = arc_is_counter_clockwise;

    tolerance = []; % Use default

    flag_arc_is_first = 0;

    [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));

    % Check size of results
    assert(isequal(size(revised_line_parameters),[1 6]));
    assert(isequal(size(revised_arc_parameters),[1 7]));

    % Check that results are correct
    assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
    assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
end

% %% BASIC test  - line first, arc last, line oriented, arc oriented, aligned, arc is negative
% fig_num = 2;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 -1];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 0; 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = []; % Use default
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% 
% %% BASIC test  - line first, arc last, line oriented, arc oriented, misaligned, arc is positive
% fig_num = 3;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 1.2];
% arc_radius               = 1;
% arc_angles               = [270 360]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 1]; 
% true_arc_is_counter_clockwise = 1;
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% %% BASIC test - line first, arc last, line oriented, arc oriented, misaligned, arc is negative
% fig_num = 4;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 -1.2];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 0;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% 
% %% BASIC test - line first, arc last, line misoriented, arc oriented, aligned, arc is positive
% fig_num = 11;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = -[1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = 0;
% line_s_end               = 1;
% 
% arc_center_xy            = [0 1];
% arc_radius               = 1;
% arc_angles               = [270 360]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 1]; 
% true_arc_is_counter_clockwise = 1;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = []; % Use default
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% %% BASIC test  - line first, arc last, line misoriented, arc oriented, aligned, arc is negative
% fig_num = 12;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = -[1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = 0;
% line_s_end               = 1;
% 
% arc_center_xy            = [0 -1];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 0;
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = []; % Use default
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% 
% %% BASIC test  - line first, arc last, line misoriented, arc oriented, misaligned, arc is positive
% fig_num = 13;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = -[1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = 0;
% line_s_end               = 1;
% 
% arc_center_xy            = [0 1.2];
% arc_radius               = 1;
% arc_angles               = [270 360]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 1]; 
% true_arc_is_counter_clockwise = 1;
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% %% BASIC test - line first, arc last, line misoriented, arc oriented, misaligned, arc is negative
% fig_num = 14;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = -[1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = 0;
% line_s_end               = 1;
% 
% arc_center_xy            = [0 -1.2];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 1; 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% %% BASIC test - line first, arc last, line oriented, arc misoriented, aligned, arc is positive
% fig_num = 21;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 1];
% arc_radius               = 1;
% arc_angles               = [360 270]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = -1;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 1]; 
% true_arc_is_counter_clockwise = 1; 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = []; % Use default
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% %% BASIC test  - line first, arc last, line oriented, arc misoriented, aligned, arc is negative
% fig_num = 22;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 -1];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 1; 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = []; % Use default
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% 
% 
% %% BASIC test  - line first, arc last, line oriented, arc misoriented, misaligned, arc is positive
% fig_num = 23;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 1.2];
% arc_radius               = 1;
% arc_angles               = [270 360]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 1;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 1]; 
% true_arc_is_counter_clockwise = 1;
% 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));
% 
% %% BASIC test - line first, arc last, line oriented, arc misoriented, misaligned, arc is negative
% fig_num = 24;
% figure(fig_num); clf;
% 
% line_unit_tangent_vector = [1 0];
% line_base_point_xy       = [0 0];
% line_s_start             = -1;
% line_s_end               = 0;
% 
% arc_center_xy            = [0 -1.2];
% arc_radius               = 1;
% arc_angles               = [90 0]*pi/180;
% arc_is_circle            = 0;
% arc_is_counter_clockwise = 0;
% 
% true_line_unit_tangent_vector = [1 0];
% true_start_point_xy = [-1 0];
% true_arc_center_xy  = [0 -1];
% true_arc_is_counter_clockwise = 1; 
% 
% % Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% line_parameters(1,1:2) = line_unit_tangent_vector;
% line_parameters(1,3:4) = line_base_point_xy;
% line_parameters(1,5)   = line_s_start;
% line_parameters(1,6)   = line_s_end;
% 
% % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
% arc_parameters(1,1:2) = arc_center_xy;
% arc_parameters(1,3)   = arc_radius;
% arc_parameters(1,4:5) = arc_angles;
% arc_parameters(1,7)   = arc_is_counter_clockwise;
% 
% tolerance = 0.4; 
% 
% flag_arc_is_first = 0;
% 
% [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
% 
% % Check size of results
% assert(isequal(size(revised_line_parameters),[1 6]));
% assert(isequal(size(revised_arc_parameters),[1 7]));
% 
% % Check that results are correct
% assert(isequal(round(revised_line_parameters,4),round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4)));
% assert(isequal(round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4), round([true_arc_center_xy arc_radius mod(arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4)));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_joinLineToArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end