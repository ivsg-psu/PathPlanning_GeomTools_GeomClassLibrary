%% script_test_fcn_geometry_fitSequentialArcs
% Exercises the function: fcn_geometry_fitSequentialArcs

% 2024_04_14 - S. Brennan
% -- wrote the code

close all;

%% Use fillArcSequence to create some test data
fig_num = 1;
figure(fig_num);
clf;

rng(1); % Fix the random number, for debugging

% arc_pattern has [1/R and L] for each segment as a row
% arc_pattern = [...
%     1/20, 15; 
%     0 20;
%     -1/5 10; 
%     0 10;
%     1/15 40; 
%     0 15
%     -1/10 20];

arc_pattern = [...
    1/20, 15; 
    0 20];

M = 10; % How many points per meter
sigma = 0.02; % The standard deviation in the points relative to the perfect function fit, in meters

[test_points, ~, ~, trueArcStartIndicies, trueNamedCurveTypes, trueParameters] = fcn_geometry_fillArcSequenceTestPoints(arc_pattern, M, sigma, -1);

% Add more noise?
if 1==0
    % Corrupt the results
    probability_of_corruption = 1;
    magnitude_of_corruption = 0.03;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end

% Add outliers?
if 1==0
    % Corrupt the results
    probability_of_corruption = 0.1;
    magnitude_of_corruption = 1;

    test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
        (probability_of_corruption), (magnitude_of_corruption), (-1));
end


% Initialize the subplots
fig_num_array(1) = fig_num;
fig_num_array(2) = 0;
fig_num_array(3) = 0;
fig_num_array(4) = 0;
figure(fig_num); clf;

% Add starter points (truth) onto subplot 2,1,1
subplot(2,2,1);
hold on;
grid on;
axis equal;
xlabel('X [meters]');
ylabel('Y [meters]');

% Plot the groups of true points
modifiedArcStartIndicies = [trueArcStartIndicies; length(test_points(:,1))];
for ith_plot = 1:length(trueArcStartIndicies(:,1))
    if ~isempty(trueNamedCurveTypes)
        current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,trueNamedCurveTypes{ith_plot},-1);
    else
        current_color = [0 0 0];
    end
    index_range = modifiedArcStartIndicies(ith_plot):modifiedArcStartIndicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end


% Perform the fit forwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 0;
[fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
    fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));


% Perform the fit backwards
fitting_tolerance = 0.1; % Units are meters
flag_fit_backwards = 1;
[fitSequence_points_backward, fitSequence_shapes_backward, fitSequence_endIndicies_backward, fitSequence_parameters_backward, fitSequence_bestFitType_backward] = ...
    fcn_geometry_fitSequentialArcs(test_points, fitting_tolerance, flag_fit_backwards, fig_num_array);

% Plot the true results
subplot(2,2,4);
fcn_geometry_plotFitSequences(trueNamedCurveTypes, trueParameters,(fig_num_array(1)));


% Compare lengths and parameters
NfitsInSequence = length(fitSequence_points_forward);

% First, make absolutely sure that the number of fits found in the forward
% direction match the same number of fits in the backward direction
if length(fitSequence_points_backward)~=NfitsInSequence
    warning('on','backtrace');
    warning('An error will be thrown at this code location as the fits were directionally different.');
    error('Found different numbers of fits in the fit sequence when comparing forward/backward directions. Code is not able to handle this yet!');
end

for ith_fit = 1:NfitsInSequence
    if ~strcmp(fitSequence_bestFitType_forward{ith_fit},fitSequence_bestFitType_backward{ith_fit})
        warning('on','backtrace');
        warning('An error will be thrown at this code location as the fits were found to be geometrically different.');
        error('Found different geometries as the best fit for the same regions when comparing forward/backward directions. Code is not able to handle this yet!');
    end
end


% Find the probable fit
fitSequence_indicies_matrix_forward = cell2mat(fitSequence_endIndicies_forward)';
fitSequence_indicies_matrix_backward = cell2mat(fitSequence_endIndicies_backward)';
probable_arc_boundary_indicies = round(mean([fitSequence_indicies_matrix_forward fitSequence_indicies_matrix_backward],2));
% probable_arc_boundary_indicies = probable_arc_boundary_indicies(1:end-1,:);

% Print and plot the results
% % Add vertical lines to indicate where the segments are TRUELY changing
% for ith_start = 1:length(arcStartIndicies)
%     plot([arcStartIndicies(ith_start) arcStartIndicies(ith_start)],[-0.1 1.1],'k-','LineWidth',5);
% end
%

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit number:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',ith_fit),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('True start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',trueArcStartIndicies(ith_fit)),10));
end
fprintf(1,'\n');

fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit start index:'),20));
for ith_fit = 1:NfitsInSequence
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.0d',probable_arc_boundary_indicies(ith_fit)),10));
end
fprintf(1,'\n');

% Add vertical lines to indicate where the segments were identified as
% changing
figure(fig_num_array(1));
subplot(2,2,3);
for ith_start = 1:NfitsInSequence
    
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_start,fitSequence_bestFitType_forward{ith_start},-1);

    plot([probable_arc_boundary_indicies(ith_start) probable_arc_boundary_indicies(ith_start)],[-0.1 1.1],'-','Color',current_color);
end

figure(fig_num_array(1));
subplot(2,2,4);

% Plot the fitted groups of points. If any of the points are mis-labeled,
% there will be one color incorrectly on top of another, for example a red
% point on top of a blue underlying point.
for ith_plot = 1:NfitsInSequence
    current_color = fcn_geometry_fillColorFromNumberOrName(ith_plot,fitSequence_bestFitType_forward{ith_plot},-1);
    index_range = probable_arc_boundary_indicies(ith_plot):probable_arc_boundary_indicies(ith_plot+1);
    plot(test_points(index_range,1),test_points(index_range,2),'.','Color',current_color,'MarkerSize',10);
end

%% Connect the fits so that the lines perfectly align with the arcs

fig_num = 23456;
figure(fig_num);clf;

revised_fitSequence_parameters_forward  = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_forward, fitSequence_parameters_forward, 0.5, fig_num);
revised_fitSequence_parameters_backward = fcn_geometry_alignGeometriesInSequence(fitSequence_bestFitType_backward,fitSequence_parameters_backward, fitting_tolerance*2, fig_num);

fcn_geometry_plotFitSequences(fitSequence_bestFitType_forward, revised_fitSequence_parameters_forward,(fig_num));
fcn_geometry_plotFitSequences(fitSequence_bestFitType_backward, revised_fitSequence_parameters_backward,(fig_num));

subplot(1,2,1);
good_axis_limits = axis;

%% Print the results


NleadCharacters = 20;
threshold = 0.15;

max_forward_error = -inf;
max_backward_error = -inf;
max_averaged_error = -inf;

fprintf(1,'\n\nPARAMETER FIT COMPARISON:\n');
for ith_fit = 1:NfitsInSequence
    fprintf(1,'\n\nFit Sequence Number: %.0d\n', ith_fit); 

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('   '),NleadCharacters));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},1)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('TRUE '),NleadCharacters));
    fcn_INTERNAL_printFitDetails(trueNamedCurveTypes{ith_fit},trueParameters{ith_fit},0)
     
    % fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD'),NleadCharacters));
    % fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit}, fitSequence_parameters_forward{ith_fit},0)
    % 
    % fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE'),NleadCharacters));
    % fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit}, fitSequence_parameters_backward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('FORWARD REV'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},revised_fitSequence_parameters_forward{ith_fit},0)

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('REVERSE REV'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_backward{ith_fit},revised_fitSequence_parameters_backward{ith_fit},0)

    averaged_parameters = (revised_fitSequence_parameters_backward{ith_fit} + revised_fitSequence_parameters_forward{ith_fit})/2;

    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('AVERAGED'),NleadCharacters));
    fcn_INTERNAL_printFitDetails(fitSequence_bestFitType_forward{ith_fit},averaged_parameters,0)
end
fprintf(1,'\n');

fprintf(1,'Max forward fitting error:  %.3f meters\n',max_forward_error);
fprintf(1,'Max backward fitting error: %.3f meters\n',max_backward_error);
fprintf(1,'Max averaged fitting error: %.3f meters\n',max_averaged_error);


%% Plot the results
comparison_fig_num = 2828;
figure(comparison_fig_num); clf;
hold on;

for ith_fit = 1:NfitsInSequence

    averaged_parameters = (revised_fitSequence_parameters_backward{ith_fit} + revised_fitSequence_parameters_forward{ith_fit})/2;

    sgtitle({'Fit quality',sprintf('Red is %.2fm error, blue is 0m error', threshold)});

    subplot(1,3,1);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_forward{ith_fit}, revised_fitSequence_parameters_forward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_forward_error = max(max_forward_error,max_error);
    title('Forward fitting');
    axis(good_axis_limits)

    subplot(1,3,2);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_backward{ith_fit}, revised_fitSequence_parameters_backward{ith_fit},...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_backward_error = max(max_backward_error,max_error);
    title('Reverse fitting');
    axis(good_axis_limits)

    subplot(1,3,3);
    curve_test_segment_length = [];    
    [flag_is_similar, points_XY_on_test_curve, minimum_distance_to_each_point, mean_error, max_error, std_dev_error] = ...
    fcn_geometry_compareCurves(...
    trueNamedCurveTypes{ith_fit}, trueParameters{ith_fit}, fitSequence_bestFitType_backward{ith_fit}, averaged_parameters,...
    (threshold), (curve_test_segment_length), (comparison_fig_num));
    max_averaged_error = max(max_averaged_error,max_error);
    title('Averaged fitting');
    axis(good_axis_limits)


end


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

%% fcn_INTERNAL_printFitDetails
function fcn_INTERNAL_printFitDetails(fit_type,fit_parameters, flag_print_header)

if contains(fit_type,{'arc','circle'})
    print_type = 'arc';            
    header_strings = {'centerX','centerY','radius','startAngle(deg)','endAngle(deg)','isCircle','turnsLeft'};
elseif contains(fit_type,{'segment','line'})
    print_type = 'line';
    header_strings = {'vect_x','vect_y','start_x','start_y','start_station','end_station'};
end

NumColumnChars = 15;

% Print header?
if flag_print_header
    % Print values
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf(' '),20));

    parameter_string = '';
    for ith_parameter = 1:length(fit_parameters)
        number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s ',header_strings{ith_parameter}),NumColumnChars+1);
        parameter_string = cat(2,parameter_string,number_string);
    end
    fprintf(1,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',parameter_string),7*NumColumnChars));
else
    % Print values
    fprintf(1,'%s',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('Fit type: %s ',print_type),20));

    % Print parameters
    parameter_string = '';
    for ith_parameter = 1:length(fit_parameters)
        if strcmp(print_type,'arc') && ((ith_parameter==4)||(ith_parameter==5))
            number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.2f ',fit_parameters(ith_parameter)*180/pi),NumColumnChars);
        else
            number_string = fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%.4f ',fit_parameters(ith_parameter)),NumColumnChars);
        end
        if strcmp(number_string(1),'-')
            parameter_string = cat(2,parameter_string,number_string, ' ');
        else
            parameter_string = cat(2,parameter_string,' ', number_string);
        end
    end
    fprintf(1,'%s\n',fcn_DebugTools_debugPrintStringToNCharacters(sprintf('%s',parameter_string),7*NumColumnChars));
end

end % Ends fcn_INTERNAL_printFitDetails









