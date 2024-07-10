%% script_test_fcn_geometry_printFitSequences
% Exercises the function: fcn_geometry_printFitSequences
% Revision history:
% 2024_04_12
% -- wrote the code

close all;

%% BASIC test - bulk line printing using defaults
flag_all_parameters_same = [];
lead_string = 'all defaults';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; %#ok<SAGROW>
    fitSequence_parameters{ith_domain}  = line_parameters; %#ok<SAGROW>

end
fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))


%% BASIC test - bulk line printing forcing header
flag_all_parameters_same = 1;
lead_string = 'forced header';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; 
    fitSequence_parameters{ith_domain}  = line_parameters; 

end
fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))

%% BASIC test - bulk line printing forcing NO header
flag_all_parameters_same = 0;
lead_string = 'forced NO header';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; 
    fitSequence_parameters{ith_domain}  = line_parameters; 

end
fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))




%% BASIC test - mixed plotting, default
flag_all_parameters_same = 0;
lead_string = 'default - mixed';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; 
    fitSequence_parameters{ith_domain}  = line_parameters; 

end

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_fitTypes{ith_domain+1} = 'arc';
fitSequence_parameters{ith_domain+1}  = arc_parameters;

fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))

%% BASIC test - mixed plotting, forced header
flag_all_parameters_same = 1;
lead_string = 'mixed - forced';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; 
    fitSequence_parameters{ith_domain}  = line_parameters; 

end

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_fitTypes{ith_domain+1} = 'arc';
fitSequence_parameters{ith_domain+1}  = arc_parameters;

fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))

%% BASIC test - mixed plotting, not forced header
flag_all_parameters_same = 0;
lead_string = 'mixed - not forced';
fid = [];

fprintf(1,'\n\n');
for ith_domain = 1:5
    line_parameters(1,1:2) = [ith_domain 0];
    line_parameters(1,3)   = ith_domain*pi/4;

    fitSequence_fitTypes{ith_domain} = 'line'; 
    fitSequence_parameters{ith_domain}  = line_parameters; 

end

% Fill the arc parameters - for listing of meaning of parameters, see fcn_geometry_fillEmptyDomainStructure
arc_center_xy            = [0 1];
arc_radius               = 1;
arc_vector_start         = [ 0 -1];
arc_vector_end           = [ 1  0];
arc_is_circle            = 0;
arc_is_counter_clockwise = 1;
arc_angles = [atan2(arc_vector_start(2),arc_vector_start(1)); atan2(arc_vector_end(2),arc_vector_end(1));];

clear arc_parameters
arc_parameters(1,1:2) = arc_center_xy;
arc_parameters(1,3)   = arc_radius;
arc_parameters(1,4:5) = arc_angles;
arc_parameters(1,6) = arc_is_circle;
arc_parameters(1,7)   = arc_is_counter_clockwise;

fitSequence_fitTypes{ith_domain+1} = 'arc';
fitSequence_parameters{ith_domain+1}  = arc_parameters;

fcn_geometry_printFitSequences(fitSequence_fitTypes, fitSequence_parameters, (flag_all_parameters_same), (lead_string), (fid))

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_printFitSequences(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end