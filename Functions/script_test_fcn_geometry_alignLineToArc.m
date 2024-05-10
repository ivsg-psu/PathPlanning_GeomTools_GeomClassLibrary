%% script_test_fcn_geometry_alignLineToArc
% Exercises the function: fcn_geometry_alignLineToArc

% Revision history:
% 2024_04_12 - Sean Brennan
% -- wrote the code
% 2024_04_19 - Sean Brennan
% -- renamed from fcn_geometry_joinLineToArc
% 2024_05_10 - Sean Brennan
% -- added test sections

close all;

%% check input orientation corrections
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
% 1

%% Basic test 1.1 - an arc nearby the line segment joined with C0 continuity
fig_num = 11;
figure(fig_num); clf;

tolerance = 0.5; % meters
shift_error = [0 0.2];

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
segment_s_start             = 0;
segment_s_end               = 1;


true_arc_center_xy  = [0 1];
% true_arc_is_counter_clockwise = 1;
% true_arc_angles     = [270 360]*pi/180;

arc_center_xy            = true_arc_center_xy+shift_error;
arc_radius               = 1;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
% arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;

continuity_level = 0;
[revised_line_parameters, revised_arc_parameters] = fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C0 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check segment results
fitted_line_params = round(revised_line_parameters,4);
true_line_params   = round([true_segment_unit_tangent_vector true_start_point_xy 0 segment_s_end-line_s_start],4);
assert(isequal(fitted_line_params, true_line_params));

% Check arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([-0.3420    0.9397    1.0000    5.0615         0         0    1.0000], fitted_arc_params));

%% Basic test 1.2 - an arc intersecting the line segment joined with C0 continuity
fig_num = 12;
figure(fig_num); clf;

tolerance = 0.5; % meters
shift_error = [-0.8 -0.2];

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
segment_s_end               = 1;


true_arc_center_xy  = [0 1];
true_arc_is_counter_clockwise = 1;
true_arc_angles     = [270 360]*pi/180;

arc_center_xy            = true_arc_center_xy+shift_error;
arc_radius               = 1;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;


continuity_level = 0;
[revised_line_parameters, revised_arc_parameters] = fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C0 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_line_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
assert(isequal([1     0    -1     0     0     1 ], fitted_line_params));

% Check the arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([-0.3420    0.9397    1.0000    5.0615         0         0    1.0000], fitted_arc_params));

%% ADVANCED test 1.9 - systematic joining line to arc, C0 continuity
fig_num = 19;

% Bit 1 (top bit)            : line precedes arc, 
% Bit 2 (2nd from top bit)   : line oriented, 
% Bit 3                      : arc oriented, 
% Bit 4                      : ends aligned, 
% Bit 5 (lowest bit)         : arc is positive
NtotalTests = (2^5);
for ith_test = 1:NtotalTests

    numDigits = 5;
    binary_string = dec2bin(ith_test-1,numDigits);

    figure(fig_num); clf;
    title_string = sprintf('Test %.0d of %.0d, %s: ',ith_test, NtotalTests, binary_string);


    % Top bit - is the line first? 0 is yes, 1 is no
    if strcmp(binary_string(1),'0')
        title_string = cat(2,title_string,'line precedes arc,');

        true_segment_unit_tangent_vector = [1 0];
        true_start_point_xy = [-1 0];
        
        segment_unit_tangent_vector = [1 0];
        segment_base_point_xy       = [-1 0];
        line_s_start             = 0;
        segment_s_end               = 1;

        arc_center_xy            = [0 1];
        arc_radius               = 1;
        arc_vector_start         = [ 0 -1];
        arc_vector_end           = [ 1  0];
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 1;

        true_arc_center_xy  = [0 1];
        true_arc_is_counter_clockwise = 1;
        true_arc_angles     = [270 360]*pi/180;

        
    else
        title_string = cat(2,title_string,'arc precedes line,');

        true_segment_unit_tangent_vector = [1 0];
        true_start_point_xy           = [0 0];

        segment_unit_tangent_vector = [1 0];
        segment_base_point_xy       = [0 0];
        line_s_start             = 0;
        segment_s_end               = 1;

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
        segment_base_point_xy = segment_base_point_xy + segment_s_end*segment_unit_tangent_vector;

        % Change the line's orientation
        segment_unit_tangent_vector = -segment_unit_tangent_vector;
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
            segment_base_point_xy = segment_base_point_xy + [0 0.2];
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
    segment_parameters(1,1:2) = segment_unit_tangent_vector;
    segment_parameters(1,3:4) = segment_base_point_xy;
    segment_parameters(1,5)   = line_s_start;
    segment_parameters(1,6)   = segment_s_end;

    % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
    arc_parameters(1,1:2) = arc_center_xy;
    arc_parameters(1,3)   = arc_radius;
    arc_parameters(1,4:5) = arc_angles;
    arc_parameters(1,7)   = arc_is_counter_clockwise;

    tolerance = 0.4; % Use default


    continuity_level = 0;
    [revised_line_parameters, revised_arc_parameters] = ...
        fcn_geometry_alignLineToArc(arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));
    sgtitle(title_string);
    pause(0.01);

    % Check size of results
    assert(isequal(size(revised_line_parameters),[1 6]));
    assert(isequal(size(revised_arc_parameters),[1 7]));

    % Check that line results
    fitted_line_params = round(revised_line_parameters,4);
    true_line_params   = round([true_segment_unit_tangent_vector true_start_point_xy 0 segment_s_end-line_s_start],4);
    assert(isequal(fitted_line_params, true_line_params));

    % Check the arc results
    fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
    true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
    assert(isequal(true_arc_params, fitted_arc_params));
end



%% Basic test 2.1 - an arc nearby the line segment joined with C1 continuity
fig_num = 21;
figure(fig_num); clf;

tolerance = 0.4; % meters
shift_error = [0 0.2];

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
segment_s_end               = 1;


true_arc_center_xy  = [0 1];
true_arc_is_counter_clockwise = 1;
true_arc_angles     = [270 360]*pi/180;

arc_center_xy            = true_arc_center_xy+shift_error;
arc_radius               = 1;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;


continuity_level = 1;
[revised_line_parameters, revised_arc_parameters] = fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C1 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_line_parameters,4);
true_line_params   = round([true_segment_unit_tangent_vector true_start_point_xy 0 segment_s_end-line_s_start],4);
assert(isequal(fitted_line_params, true_line_params));

% Check the arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal(true_arc_params, fitted_arc_params));

%% Basic test 2.2 - an arc intersecting the line segment joined with C1 continuity
fig_num = 22;
figure(fig_num); clf;

tolerance = 1; % meters
shift_error = [-0.8 -0.2];

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
segment_s_end               = 1;


true_arc_center_xy  = [0 1];
true_arc_is_counter_clockwise = 1;
true_arc_angles     = [270 360]*pi/180;

arc_center_xy            = true_arc_center_xy+shift_error;
arc_radius               = 1;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;


continuity_level = 1;
[revised_line_parameters, revised_arc_parameters] = fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C0 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_line_parameters,4);
true_line_params   = round([true_segment_unit_tangent_vector true_start_point_xy 0 segment_s_end-line_s_start],4);
assert(isequal(round(true_line_params,4), fitted_line_params));

% Check the arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal(round(true_arc_params,4), fitted_arc_params));


%% ADVANCED test 2.9 - joining line to arc, C1 continuity
fig_num = 29;

% Bit 1 (top bit)            : line precedes arc, 
% Bit 2 (2nd from top bit)   : line oriented, 
% Bit 3                      : arc oriented, 
% Bit 4                      : ends aligned, 
% Bit 5 (lowest bit)         : arc is positive
NtotalTests = (2^5);
for ith_test = 1:NtotalTests

    numDigits = 5;
    binary_string = dec2bin(ith_test-1,numDigits);

    figure(fig_num); clf;
    title_string = sprintf('Test %.0d of %.0d, %s: ',ith_test, NtotalTests, binary_string);
    title(title_string)

    % Top bit - is the line first? 0 is yes, 1 is no
    if strcmp(binary_string(1),'0')
        title_string = cat(2,title_string,'line precedes arc,');

        true_segment_unit_tangent_vector = [1 0];
        true_start_point_xy = [-1 0];
        
        segment_unit_tangent_vector = [1 0];
        segment_base_point_xy       = [-1 0];
        line_s_start             = 0;
        segment_s_end               = 1;

        arc_center_xy            = [0 1];
        arc_radius               = 1;
        arc_vector_start         = [ 0 -1];
        arc_vector_end           = [ 1  0];
        arc_is_circle            = 0;
        arc_is_counter_clockwise = 1;

        true_arc_center_xy  = [0 1];
        true_arc_is_counter_clockwise = 1;
        true_arc_angles     = [270 360]*pi/180;

        
    else
        title_string = cat(2,title_string,'arc precedes line,');

        true_segment_unit_tangent_vector = [1 0];
        true_start_point_xy           = [0 0];

        segment_unit_tangent_vector = [1 0];
        segment_base_point_xy       = [0 0];
        line_s_start             = 0;
        segment_s_end               = 1;

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
        segment_base_point_xy = segment_base_point_xy + segment_s_end*segment_unit_tangent_vector;

        % Change the line's orientation
        segment_unit_tangent_vector = -segment_unit_tangent_vector;
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
            segment_base_point_xy = segment_base_point_xy + [0 0.2];
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
    segment_parameters(1,1:2) = segment_unit_tangent_vector;
    segment_parameters(1,3:4) = segment_base_point_xy;
    segment_parameters(1,5)   = line_s_start;
    segment_parameters(1,6)   = segment_s_end;

    % Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
    arc_parameters(1,1:2) = arc_center_xy;
    arc_parameters(1,3)   = arc_radius;
    arc_parameters(1,4:5) = arc_angles;
    arc_parameters(1,7)   = arc_is_counter_clockwise;

    tolerance = 0.4; % Use default


    continuity_level = 1;
    [revised_line_parameters, revised_arc_parameters] = ...
        fcn_geometry_alignLineToArc(arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));
    sgtitle(title_string);
    pause(0.01);

    % Check size of results
    assert(isequal(size(revised_line_parameters),[1 6]));
    assert(isequal(size(revised_arc_parameters),[1 7]));

    % Check that line results
    fitted_line_params = round(revised_line_parameters,4);
    true_line_params   = round([true_segment_unit_tangent_vector true_start_point_xy 0 segment_s_end-line_s_start],4);
    assert(isequal(fitted_line_params, true_line_params));

    % Check the arc results
    fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
    true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
    assert(isequal(true_arc_params, fitted_arc_params));
end

%% Test that threshold will not cause fitting if error is too large
fig_num = 2;
figure(fig_num); clf;

tolerance = 0.01; % Use very low tolerance to force fit to NOT occur
shift_error = [0 0.2];

title_string = 'Checking that a large error causes fit to not occur';

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = [1 0];
segment_base_point_xy       = [-1 0] + shift_error;
line_s_start             = 0;
segment_s_end               = 1;

arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


true_arc_center_xy  = [0 1];
true_arc_is_counter_clockwise = 1;
true_arc_angles     = [270 360]*pi/180;

% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;


continuity_level = 1;
[revised_line_parameters, revised_arc_parameters] = fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

% Check size of results
assert(isempty(revised_line_parameters));
assert(isempty(revised_arc_parameters));


%% Basic test 3.1 - an arc nearby the line segment joined with C2 continuity
fig_num = 31;
figure(fig_num); clf;

tolerance = 0.4; % meters
shift_error = [0 0.1];

true_segment_unit_tangent_vector = [1 0];
true_start_point_xy = [-1 0];

segment_unit_tangent_vector = true_segment_unit_tangent_vector;
segment_base_point_xy       = true_start_point_xy;
line_s_start             = 0;
segment_s_end               = 1;


true_arc_center_xy  = [0 1];
true_arc_is_counter_clockwise = 1;
true_arc_angles     = [270 360]*pi/180;

arc_center_xy            = true_arc_center_xy+shift_error;
arc_radius               = 1;
arc_vector_start         = [cos(-70*pi/180) sin(-70*pi/180)];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];


% Get the line fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
segment_parameters(1,1:2) = segment_unit_tangent_vector;
segment_parameters(1,3:4) = segment_base_point_xy;
segment_parameters(1,5)   = line_s_start;
segment_parameters(1,6)   = segment_s_end;

% Get the arc fit details from parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,7)   = arc_is_counter_clockwise;


continuity_level = 2;
[revised_line_parameters, revised_arc_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C2 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_line_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));

% Check the arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));

% Check the spiral results
assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));

%% Checking case that fails

fig_num = 99;
figure(fig_num); clf;


segment_parameters = [ 0.7317    0.6816    3.6549   -3.9239   15.7332   33.6332];
arc_parameters =  [-0.1225   20.3593   20.3705    4.7184    5.5674         0    1.0000];

flag_arc_is_first = 1;
tolerance = 0.5;         % meters
continuity_level = 1;
[revised_line_parameters, revised_arc_parameters, revised_spiral_join_parameters] = ...
    fcn_geometry_alignLineToArc(...
    arc_parameters, segment_parameters, (tolerance), (continuity_level), (fig_num));

title('Checking that arc is joined to the line: C1 continuous');


% Check size of results
assert(isequal(size(revised_line_parameters),[1 6]));
assert(isequal(size(revised_arc_parameters),[1 7]));

% Check that line results
fitted_line_params = round(revised_line_parameters,4);
% true_line_params   = round([true_line_unit_tangent_vector true_start_point_xy 0 line_s_end-line_s_start],4);
assert(isequal([1.0000         0   -1.0000         0         0    0.2326], round(fitted_line_params,4)));

% Check the arc results
fitted_arc_params = round([revised_arc_parameters(1,1:3) mod(revised_arc_parameters(1,4:5),2*pi) revised_arc_parameters(1,6:7)],4);
% true_arc_params = round([true_arc_center_xy arc_radius mod(true_arc_angles,2*pi) arc_is_circle true_arc_is_counter_clockwise],4);
assert(isequal([ 0    1.1000    1.0000    5.4955         0         0    1.0000], round(fitted_arc_params,4)));

% Check the spiral results
assert(isequal([1.5662         0   -0.7674         0         0    1.0000], round(revised_spiral_join_parameters,4)));


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_alignLineToArc(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end