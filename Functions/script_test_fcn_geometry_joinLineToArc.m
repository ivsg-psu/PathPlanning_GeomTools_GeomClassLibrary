%% script_test_fcn_geometry_joinLineToArc
% Exercises the function: fcn_geometry_joinLineToArc
% Revision history:
% 2024_04_12
% -- wrote the code

close all;

%% BASIC test

% Bit 1 (top bit)            : line precedes arc, 
% Bit 2 (2nd from top bit)   : line oriented, 
% Bit 3                      : arc oriented, 
% Bit 4                      : ends aligned, 
% Bit 5 (lowest bit)         : arc is positive
NtotalTests = (2^5);
for ith_test = 1:NtotalTests

    numDigits = 5;
    binary_string = dec2bin(ith_test-1,numDigits);

    fig_num = 1;
    figure(fig_num); clf;
    title_string = sprintf('Test %.0d of %.0d, %s: ',ith_test, NtotalTests, binary_string);


    % Top bit - is the line first? 0 is yes, 1 is no
    if strcmp(binary_string(1),'0')
        title_string = cat(2,title_string,'line precedes arc,');

        true_line_unit_tangent_vector = [1 0];
        true_start_point_xy = [-1 0];
        
        line_unit_tangent_vector = [1 0];
        line_base_point_xy       = [-1 0];
        line_s_start             = 0;
        line_s_end               = 1;

        arc_center_xy            = [0 1];
        arc_radius               = 1;
        arc_vector_start         = [ 0 -1];
        arc_vector_end           = [ 1  0];
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 1;

        true_arc_center_xy  = [0 1];
        true_arc_is_counter_clockwise = 1;
        true_arc_angles     = [270 360]*pi/180;

        flag_arc_is_first = 0;
    else
        title_string = cat(2,title_string,'arc precedes line,');

        true_line_unit_tangent_vector = [1 0];
        true_start_point_xy           = [0 0];

        line_unit_tangent_vector = [1 0];
        line_base_point_xy       = [0 0];
        line_s_start             = 0;
        line_s_end               = 1;

        arc_center_xy            = [0 1];
        arc_radius               = 1;
        arc_vector_start         = [-1  0];
        arc_vector_end           = [ 0 -1];  
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 1;

        true_arc_center_xy  = [0 1];
        true_arc_is_counter_clockwise = 1;
        true_arc_angles          = [180 270]*pi/180;
        flag_arc_is_first = 1;
    end

    % Next from top bit - is the line oriented correctly? 0 is yes, 1 is no
    if strcmp(binary_string(2),'0')
        title_string = cat(2,title_string,'line oriented,');
    else
        title_string = cat(2,title_string,'line misoriented,');
        % Move the start point of the line to the end
        line_base_point_xy = line_base_point_xy + line_s_end*line_unit_tangent_vector;

        % Change the line's orientation
        line_unit_tangent_vector = -line_unit_tangent_vector;
    end

    % Next from top bit - is the arc oriented correctly? 0 is yes, 1 is no
    if strcmp(binary_string(3),'0')
        title_string = cat(2,title_string,'arc oriented,');
    else
        title_string = cat(2,title_string,'arc misoriented,');
        temp = arc_vector_start;
        arc_vector_start = arc_vector_end;
        arc_vector_end   = temp;
        arc_is_counter_clockwise = ~arc_is_counter_clockwise;
    end

    % Next from top bit - are the ends aligned? 0 is yes, 1 is no
    if strcmp(binary_string(4),'0')
        title_string = cat(2,title_string,'ends aligned,');
    else
        title_string = cat(2,title_string,'ends misaligned,');
        if strcmp(binary_string(1),'0')
            % Line is first, misalign the arc
            arc_center_xy = arc_center_xy+[0 0.2];
        else
            % Arc is first, misalign the line
            line_base_point_xy = line_base_point_xy + [0 0.2];
        end
    end

    % Bottom bit - is the arc oriented to the left? 0 is yes, 1 is no
    if strcmp(binary_string(5),'0')
        title_string = cat(2,title_string,'arc left,');
    else
        title_string = cat(2,title_string,'arc right,');

        % Flip the y values
        arc_vector_start(2) = -arc_vector_start(2);
        arc_vector_end(2)   = -arc_vector_end(2);
        arc_is_counter_clockwise = ~arc_is_counter_clockwise;
        arc_center_xy = -1*arc_center_xy;
        
        true_arc_center_xy = -1*true_arc_center_xy;        
        true_arc_is_counter_clockwise = 0;
        true_arc_angles = 2*pi-true_arc_angles;
    end

    arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

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

    tolerance = 0.4; % Use default



    [revised_line_parameters, revised_arc_parameters] = fcn_geometry_joinLineToArc(line_parameters, arc_parameters, flag_arc_is_first, (tolerance),(fig_num));
    sgtitle(title_string);
    pause(0.01);

    % Check size of results
    assert(isequal(size(revised_line_parameters),[1 6]));
    assert(isequal(size(revised_arc_parameters),[1 7]));

    % Check that line results
    fitted_line_params = round(revised_line_parameters,4);
    true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
    assert(isequal(fitted_line_params, true_line_params));

    % Check the arc results
    fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
    true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
    assert(isequal(true_arc_params, fitted_arc_params));
end


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_joinLineToArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end