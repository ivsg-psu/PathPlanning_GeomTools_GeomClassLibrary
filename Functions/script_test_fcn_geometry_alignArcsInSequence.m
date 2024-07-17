%% script_test_fcn_geometry_alignArcsInSequence
% Tests the function fcn_geometry_alignArcsInSequence

% Revision history:
% 2024_04_19 - S Brennan
% -- wrote the code
% 2024_05_02 - S Brennan / A. Batchu
% -- add more test cases
% 2024_06_26 - S Brennan
% -- working with test track data

close all;

%% REAL WORLD TEST CASE
% Test with real-world data (test track)

flag_force_load = 1;

% Check to see if the data was loaded earlier
if 1==flag_force_load || ~exist('XY_data','var')

    % Check to see if XY data for the centerline of the original track lane was loaded earlier
    mat_filename = fullfile(cd,'Data','Centerline_OriginalTrackLane_InnerMarkerClusterCenterOfDoubleYellow.mat');
    if exist(mat_filename,'file')
        load(mat_filename,'XY_data');
    end
end

Nfittings = 0;
clear fit_results

% Fitting parameters are:
% flag_fit_type S-gap  t-gap  keep_every 
fitting_parameters_to_test = [...
    0 30 2.0 100
    1 30 2.0 100
    ];


for ith_fit = 1:length(fitting_parameters_to_test(:,1))

    flag_fit_type = fitting_parameters_to_test(ith_fit,1);
    S_tolerance   = fitting_parameters_to_test(ith_fit,2);
    t_tolerance   = fitting_parameters_to_test(ith_fit,3);
    keep_every    = fitting_parameters_to_test(ith_fit,4);

    % Check to see if the fits were calculated earlier
    if 1==flag_force_load || ~exist('fitSequence_points_forward','var') || ~exist('fitSequence_points_backward','var')

        fitting_tolerance = [S_tolerance t_tolerance]; % Units are meters, in St form

        % Since the XY data is very dense, keep only some of points
        indicies = (1:length(XY_data(:,1)))';
        small_XY_data_indicies = find(0==mod(indicies,keep_every));
        small_XY_data = XY_data(small_XY_data_indicies,:);


        %% Fit the data - one can choose many options
        % Option 0 is to fit forward sequentially, 
        % Option 1 is to fit backward sequentially,
        % Option 2 is to fit using Hough (only works for small data),
        % Option 3 is to fit using curvature,
        if 0==flag_fit_type || 1==flag_fit_type
            %%%%%
            % Perform the fit forwards or backwards
            fig_num = ith_fit;
            figure(fig_num);
            clf;

            flag_fit_backwards = flag_fit_type;
            [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters, fitSequence_bestFitType] = ...
                fcn_geometry_fitSequentialArcsC2(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);
            title('Forwards sequential fit');
        end

        if 2==flag_fit_type
            % %%%%%
            % % HOUGH FIT
            % points = small_XY_data;
            % transverse_tolerance = fitting_tolerance(1,2);
            % station_tolerance    = fitting_tolerance(1,1);
            % points_required_for_agreement = 4;
            % flag_force_circle_fit = 0;
            % expected_radii_range = [];
            % flag_find_only_best_agreement = 0;
            % flag_use_permutations = [];
            % % Perform the fit using Hough
            % fig_num = 11;
            % figure(fig_num);
            % clf;
            %
            % domains  = ...
            % fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
            %         (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));
            % title('Hough fit');
            %
            % % [fitSequence_points_forward, fitSequence_shapes_forward, fitSequence_endIndicies_forward, fitSequence_parameters_forward, fitSequence_bestFitType_forward] = ...
            % %     fcn_geometry_fitSequentialArcs(small_XY_data, fitting_tolerance, flag_fit_backwards, fig_num);

        end

        

    end

    %% Perform alignment of forward arcs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %           _ _
    %     /\   | (_)                  /\
    %    /  \  | |_  __ _ _ __       /  \   _ __ ___ ___
    %   / /\ \ | | |/ _` | '_ \     / /\ \ | '__/ __/ __|
    %  / ____ \| | | (_| | | | |   / ____ \| | | (__\__ \
    % /_/    \_\_|_|\__, |_| |_|  /_/    \_\_|  \___|___/
    %                __/ |
    %               |___/
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Align%20%20Arcs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Connect the fits so that the lines perfectly align with the arcs
    fig_num = 100*ith_fit;
    figure(fig_num);
    clf;
    continuity_level = 2;

    fits_to_check_types = fitSequence_bestFitType;
    fits_to_check_parameters = fitSequence_parameters;

    [revised_fitSequence_types, revised_fitSequence_parameters, max_feasibility_distance, flag_fit_failed] =  ...
        fcn_geometry_alignArcsInSequence(fits_to_check_types, fits_to_check_parameters, fitting_tolerance, (continuity_level), (fig_num));

    if 0==flag_fit_failed
        data_to_analyze = small_XY_data; % small_XY_data
        plot(data_to_analyze(:,1),data_to_analyze(:,2),'k.','MarkerSize',5);

        % Check the fits
        fig_num = 1000*ith_fit;
        figure(fig_num);
        clf;

        threshold           = []; %max(max_feasibility_distance,fitting_tolerance(1,2));
        curve_test_segment_length = 0.5; % Check every 0.5 meters;

        % [flag_is_similar, minimum_distance_to_each_point, mean_error_forward, max_error_forward, std_dev_error_forward] = ...
        [flag_is_similar, minimum_distance_to_each_point, indicies_of_nearest_reference_points, mean_error, max_error, std_dev_error] = ...
            fcn_geometry_comparePointsToCurve(...
            revised_fitSequence_types, revised_fitSequence_parameters, data_to_analyze, ...
            (threshold), (curve_test_segment_length), (fig_num));
    else

        revised_fitSequence_types = {};
        revised_fitSequence_parameters = {};
        mean_error = inf;
        max_error  = inf;
        std_dev_error = inf;
    end

    %% Save results
    Nfittings = Nfittings+1;

    fit_results{Nfittings}.flag_fit_type                    = flag_fit_type; %#ok<SAGROW>
    fit_results{Nfittings}.fit_failed                       = flag_fit_failed; %#ok<SAGROW>
    fit_results{Nfittings}.S_tolerance                      = S_tolerance; %#ok<SAGROW>
    fit_results{Nfittings}.t_tolerance                      = t_tolerance; %#ok<SAGROW>
    fit_results{Nfittings}.keep_every                       = keep_every; %#ok<SAGROW>
    fit_results{Nfittings}.fitSequence_bestFitType          = fitSequence_bestFitType; %#ok<SAGROW>
    fit_results{Nfittings}.fitSequence_parameters           = fitSequence_parameters; %#ok<SAGROW>
    fit_results{Nfittings}.fitSequence_bestFitType          = fitSequence_bestFitType; %#ok<SAGROW>
    fit_results{Nfittings}.revised_fitSequence_types        = revised_fitSequence_types; %#ok<SAGROW>
    fit_results{Nfittings}.revised_fitSequence_parameters   = revised_fitSequence_parameters; %#ok<SAGROW>
    fit_results{Nfittings}.mean_error                       = mean_error; %#ok<SAGROW>
    fit_results{Nfittings}.max_error                        = max_error; %#ok<SAGROW>
    fit_results{Nfittings}.std_dev_error                    = std_dev_error; %#ok<SAGROW>

    %% Summarize the fit
    fprintf(1,'\n');
    fprintf(1,'SUMMARY for fit number: %.0f\n',Nfittings);
    fprintf(1,'Fit type? %.0f\n',fit_results{Nfittings}.flag_fit_type);   
    fprintf(1,'Failed forward? %.0f\n',fit_results{Nfittings}.fit_failed);
    fprintf(1,'Station tolerance: %.3fm\n',fit_results{Nfittings}.S_tolerance);
    fprintf(1,'Transverse tolerance: %.3fm\n',fit_results{Nfittings}.t_tolerance);
    fprintf(1,'Keep every # points: %.0f\n',fit_results{Nfittings}.keep_every);
    fprintf(1,'\n');
    fprintf(1,'Number of forward fits: %.0d\n',length(fit_results{Nfittings}.revised_fitSequence_types));
    fprintf(1,'Mean forward fitting error:  %.3f meters\n',fit_results{Nfittings}.mean_error);
    fprintf(1,'Max forward fitting error:  %.3f meters\n',fit_results{Nfittings}.max_error);
    fprintf(1,'Std forward fitting error:  %.3f meters\n', fit_results{Nfittings}.std_dev_error);
    fprintf(1,'\n');
    fprintf(1,'ORIGINAL fits:\n');
    fcn_geometry_printFitSequences(fit_results{Nfittings}.fitSequence_bestFitType, fit_results{Nfittings}.fitSequence_parameters , (1), ('Original fit'), (1))
    fprintf(1,'ALIGNED fits:\n');
    fcn_geometry_printFitSequences(fit_results{Nfittings}.revised_fitSequence_types, fit_results{Nfittings}.revised_fitSequence_parameters , (1), ('Aligned fit'), (1))


end