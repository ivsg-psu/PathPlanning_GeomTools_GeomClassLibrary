%% script_test_fcn_geometry_estimateSpiralLength
% Tests the function: fcn_geometry_estimateSpiralLength

% Revision history:
% 2024_07_30 - S. Brennan
% -- wrote the code

close all

%% Simple fit test
fig_num = 1;
figure(fig_num);
clf;

circle1_parameters = [0 1  1];
circle2_radius = 0.99;
offset = 0.001;
circle2_parameters = [0 circle2_radius+offset  circle2_radius];

[estimated_spiralLength, estimated_spiralStartAngle, estimated_spiralEndAngle, angle_larger_to_smaller] = ...
    fcn_geometry_estimateSpiralLength(circle1_parameters, circle2_parameters, fig_num);

% Calculate the true values
spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], -1);

h0           = spiral_join_parameters(1,3);
spiralLength = spiral_join_parameters(1,4);
% x0           = spiral_join_parameters(1,3);
% y0           = spiral_join_parameters(1,4);
K0           = spiral_join_parameters(1,5);
Kf           = spiral_join_parameters(1,6);
analytical_end_angle   = h0 + (Kf+K0)/2*spiralLength;

percent_agreement = 0.05;

% Check that length is within percent agreement
assert(spiralLength/estimated_spiralLength<(1+percent_agreement))
assert(spiralLength/estimated_spiralLength>(1-percent_agreement))

% Check that startAngle is within percent agreement
assert(h0/estimated_spiralStartAngle<(1+percent_agreement))
assert(h0/estimated_spiralStartAngle>(1-percent_agreement))

% Check that length is within percent agreement
assert(analytical_end_angle/estimated_spiralEndAngle<(1+percent_agreement))
assert(analytical_end_angle/estimated_spiralEndAngle>(1-percent_agreement))

%% Compare times between estimation versus actual calculation

circle1_parameters = [0 1  1];
circle2_radius = 0.99;
offset = 0.001;
circle2_parameters = [0 circle2_radius+offset  circle2_radius];


%%%%%%


% Perform the calculation in slow mode
% fig_num = [];
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], -1);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [estimated_spiralLength, estimated_spiralStartAngle, estimated_spiralEndAngle, angle_larger_to_smaller] = ...
        fcn_geometry_estimateSpiralLength(circle1_parameters, circle2_parameters, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_findAgreementsOfPointsToCircle:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

% Find the correct values
h0           = spiral_join_parameters(1,3);
spiralLength = spiral_join_parameters(1,4);
% x0           = spiral_join_parameters(1,3);
% y0           = spiral_join_parameters(1,4);
K0           = spiral_join_parameters(1,5);
Kf           = spiral_join_parameters(1,6);
analytical_end_angle   = h0 + (Kf+K0)/2*spiralLength;


percent_agreement = 0.05;

% Check that length is within percent agreement
assert(spiralLength/estimated_spiralLength<(1+percent_agreement))
assert(spiralLength/estimated_spiralLength>(1-percent_agreement))

% Check that startAngle is within percent agreement
assert(h0/estimated_spiralStartAngle<(1+percent_agreement))
assert(h0/estimated_spiralStartAngle>(1-percent_agreement))

% Check that length is within percent agreement
assert(analytical_end_angle/estimated_spiralEndAngle<(1+percent_agreement))
assert(analytical_end_angle/estimated_spiralEndAngle>(1-percent_agreement))

%% Calculate fits
flag_do_plots = 0;

% Loop through offsets
Noffsets = 20;
Ncircles = 30;
offsets_to_try = flipud(logspace(log10(0.01/1000),log10(0.1),Noffsets)');

fraction_of_radial_space_for_r2 = 1 - logspace(log10(0.01/1000),log10(0.2),Ncircles)';
fraction_of_radial_space_for_r2 = flipud(fraction_of_radial_space_for_r2);

% Initialize arrays that save results
length_slopes = nan*offsets_to_try;
length_intercepts = nan*offsets_to_try;
startAngle_slopes = nan*offsets_to_try;
startAngle_intercepts = nan*offsets_to_try;
clear fit_points_x fit_points_y
fit_points_y{Noffsets} = [];
fit_points_x{Noffsets} = [];


for ith_offset = 1:length(offsets_to_try)
    fprintf(1,'Testing offset: %.0d of %.0d\n',ith_offset,Noffsets)
    offset = offsets_to_try(ith_offset,1);

    % The amount of space within the radius of the unit-radius circle must contain
    % the offset gap and the remainder maximum size of the radius. Usually,
    % the fraction is close to 1, just slightly less.
    circle2_radii = (ones(Ncircles,1)-offset).*(fraction_of_radial_space_for_r2);

    lengths = nan*circle2_radii;
    start_angles = nan*circle2_radii;
    end_angles = nan*circle2_radii;

    if 1==flag_do_plots
        fig_num = 111111;
        figure(fig_num); clf;
    end

    for jth_radius = 1:Ncircles

        fprintf(1,'\tTesting radius: %.0d of %.0d, ',jth_radius,Ncircles);
        circle2_radius = circle2_radii(jth_radius,1);

        % If the difference is on the order of the numerical tolerance of
        % MATLAB, stop the calculations
        if log(1-circle2_radius)<-15
            break
        end

        circle1_parameters = [0 1  1];
        circle2_parameters = [0 circle2_radius+offset  circle2_radius];

        % Plot the circles
        if 1==flag_do_plots

            figure(383838);
            clf; hold on;
            axis equal
            grid on;
            fcn_geometry_plotGeometry('circle',circle1_parameters,0.01);
            fcn_geometry_plotGeometry('circle',circle2_parameters,0.01);
            title(sprintf('Testing radius: %.0d of %.0d',jth_radius,Ncircles));
            pause(0.02);
        end

        if 1 == 1
            if 1==flag_do_plots
                spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], fig_num);
            else
                spiral_join_parameters = fcn_geometry_spiralFromCircleToCircle(circle1_parameters, circle2_parameters, [], -1);
            end
            h0           = spiral_join_parameters(1,3);
            spiralLength = spiral_join_parameters(1,4);
            x0           = spiral_join_parameters(1,3);
            y0           = spiral_join_parameters(1,4);
            K0           = spiral_join_parameters(1,5);
            Kf           = spiral_join_parameters(1,6);
            analytical_end_angle   = h0 + (Kf+K0)/2*spiralLength;

            % Ignore any solutions that wrap around more than 180. We never see
            % these on roads
            if spiralLength>=pi
                break
            end

            % If the length is on the order of the numerical tolerance of
            % MATLAB, remove this length from the calculations
            if log(spiralLength)<-10
                spiralLength = nan;
            end


            lengths(jth_radius,1) = spiralLength;
            start_angles(jth_radius,1) = -h0;
            end_angles(jth_radius,1) = analytical_end_angle;
            
            % Compare the actual and predicted length
            % [estimated_spiralLength, estimated_spiralStartAngle, estimated_spiralEndAngle, angle_larger_to_smaller] = ...
            %    fcn_geometry_estimateSpiralLength(circle1_parameters, circle2_parameters, fig_num);

            estimated_length = fcn_INTERNAL_estimateSpiralLength(circle2_radius, offset);
            estimated_startAngle = fcn_INTERNAL_estimateSpiralStartAngle(circle2_radius, offset);


            fprintf(1,'\tOffset: %.6f, Radius: %.6f, Length Actual: %.6f  Fit: %.6f, startAngle Actual: %.6f  Fit: %.6f,  lenRatio: %.4f angRatio: %.4f \n', ...
                offset, circle2_radius, ...
                spiralLength, estimated_length, ...
                -h0, estimated_startAngle, ...
                spiralLength/estimated_length, -h0/estimated_startAngle);
        end
    end

    % Remove NaN values
    good_indicies = ~isnan(lengths);
    circle2_radii = circle2_radii(good_indicies,:);
    start_angles = start_angles(good_indicies,:);
    end_angles = end_angles(good_indicies,:);
    lengths = lengths(good_indicies,:);

    if 1==flag_do_plots

        figure(22020);
        clf;
        loglog(1-circle2_radii, abs(start_angles*180/pi),'r.-');
        xlabel('log(1- circle2radius)')
        ylabel('log(start angles in deg)')


        % figure(22021);
        % clf;
        % loglog(1-circle2_radii, lengths,'.-');
        % xlabel('log(1-circle2_radius)')
        % ylabel('log(length)')
    end

    %%%%%
    % Fit the length
    % Check the linear regression
    x_data = log(1-circle2_radii);
    y_data = log(lengths);

    fit_points_x{ith_offset} = x_data;
    fit_points_y{ith_offset} = y_data;


    % Fit equation of form
    % y = x*m+b --> Y = [X 1]*[m b]
    X_vec = [x_data ones(length(x_data(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*y_data);
    length_slopes(ith_offset,1) = fit_params(1,1);
    length_intercepts(ith_offset,1) = fit_params(2,1);

    y_pred = [x_data ones(length(x_data(:,1)),1)]*fit_params;

    if 1==flag_do_plots

        figure(37373);
        clf;
        hold on;
        plot(x_data,y_data,'k.','MarkerSize',10);
        plot(x_data,y_pred,'b-');
        xlabel('log(1-circle2radius)')
        ylabel('log(length)')
    end

    %%%%%
    % Fit the start angle
    % Check the linear regression
    x_data = log(1-circle2_radii);
    y_data = log(start_angles);

    % Fit equation of form
    % y = x*m+b --> Y = [X 1]*[m b]
    X_vec = [x_data ones(length(x_data(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*y_data);
    startAngle_slopes(ith_offset,1) = fit_params(1,1);
    startAngle_intercepts(ith_offset,1) = fit_params(2,1);

    y_pred = [x_data ones(length(x_data(:,1)),1)]*fit_params;

    if 1==1  %flag_do_plots

        figure(66665);
        clf;
        hold on;
        plot(x_data,y_data,'k.','MarkerSize',10);
        plot(x_data,y_pred,'b-');
        xlabel('log(1-circle2radius)')
        ylabel('log(startAngles)')
    end
    
end
save('Data\SpiralApproximation','length_intercepts','length_slopes','offsets_to_try','startAngle_slopes','startAngle_intercepts');

%% Check how the fit parameters change versus inputs
load('Data\SpiralApproximation','length_intercepts','length_slopes','offsets_to_try','startAngle_slopes','startAngle_intercepts');


%%
% Fit the slope for lengths
if 1==0
    figure(47474);
    clf;
    hold on;
    % Fit the data?
    xdata = log(-log(offsets_to_try));
    ydata = log(0.5-length_slopes);

    % Fit equation of form
    % y = x2*a + x*b + c --> Y = [X^3 X^2 X 1]*[a b c]'
    X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*ydata);
    fprintf(1,'\n\nSLOPE FITTING: \n')
    fprintf(1,'Fit parameters for y = x^3 + x^2 + x + 1 for y=log(0.5-slope) and x=log(-log(offset)):\n')
    disp(fit_params);

    y_pred = X_vec*fit_params;
    plot(xdata,y_pred,'r.-','MarkerSize',10,'LineWidth',3);
    plot(xdata,ydata, 'b.','MarkerSize',20,'LineWidth',3);
    xlabel('log(log(offsets))');
    ylabel('log(0.5 - slopes)');

    fprintf(1,'Comparison of actual and fit slopes:\n');
    for ith_slope = 1:length(length_slopes)
        offset = offsets_to_try(ith_slope,1);
        fit_slope = fcn_INTERNAL_estimateLengthSlopeFromOffset(offset);
        fprintf(1,'Offset: %.9f, lengthIntercept Actual: \t %.9f  Fit: \t%.9f Ratio: %.4f\n',offset, length_slopes(ith_slope,1),fit_slope, length_slopes(ith_slope,1)/fit_slope);
    end
end

%%
% Fit the intercepts for lengths

if 1==0
    figure(656565);
    clf;
    hold on;
    grid on;

    % Fit the data?
    xdata = log(offsets_to_try);
    ydata = length_intercepts;

    % Fit equation of form
    % y = x2*a + x*b + c --> Y = [X^3 X^2 X 1]*[a b c]'
    X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];
    % X_vec = [xdata.^2 xdata ones(length(xdata(:,1)),1)];
    % X_vec = [xdata ones(length(xdata(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*ydata);
    fprintf(1,'\n\nINTERCEPT FITTING: \n')

    fprintf(1,'Fit parameters for y = x^3 + x^2 + x + 1 for y=log(0.5-slope) and x=log(-log(offset)):\n')
    disp(fit_params);

    y_pred = X_vec*fit_params;
    plot(xdata,y_pred,'r.-','MarkerSize',20,'LineWidth',3);
    plot(xdata,ydata, 'b.','MarkerSize',10,'LineWidth',3);
    xlabel('log(offsets)');
    ylabel('intercepts');

    fprintf(1,'Comparison of actual and fit intercepts:\n');
    for ith_intercept = 1:length(length_intercepts)
        offset = offsets_to_try(ith_intercept,1);
        fit_intercept = fcn_INTERNAL_estimateLengthInterceptFromOffset(offset);
        fprintf(1,'Offset: %.9f, lengthIntercept Actual: \t %.9f  Fit: \t%.9f  Ratio: %.4f\n',offset, length_intercepts(ith_intercept,1),fit_intercept, length_intercepts(ith_intercept,1)/fit_intercept);
    end
end

%% startAngleSlopes
% Fit the slope for startAngles
if 1==1
    figure(888878);
    clf;
    hold on;

    % Fit the data?
    xdata = log(-log(offsets_to_try));
    ydata = startAngle_slopes;

    % Fit equation of form
    % y = x2*a + x*b + c --> Y = [X^3 X^2 X 1]*[a b c]'
    X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*ydata);
    fprintf(1,'\n\nSLOPE FITTING for startAngles: \n')
    fprintf(1,'Fit parameters for y = x^3 + x^2 + x + 1 for y=startAngleSlopes and x=log(-log(offset)):\n')
    disp(fit_params);

    y_pred = X_vec*fit_params;
    plot(xdata,y_pred,'r.-','MarkerSize',10,'LineWidth',3);
    plot(xdata,ydata, 'b.','MarkerSize',20,'LineWidth',3);
    xlabel('log(-log(offsets))');
    ylabel('startAngle_slopes');
    
    fprintf(1,'\n\nComparison of actual and fit startAngleSlopes:\n');
    for ith_slope = 1:length(startAngle_slopes)
        offset = offsets_to_try(ith_slope,1);
        fit_slope = fcn_INTERNAL_estimateStartAngleSlopeFromOffset(offset);
        fprintf(1,'Offset: %.9f: startAngleSlopes: Actual: \t %.9f  Fit: \t%.9f Ratio: %.4f \n',offset, startAngle_slopes(ith_slope,1),fit_slope, startAngle_slopes(ith_slope,1)/fit_slope);
    end
end

%% startAngleIntercept
% Fit the intercepts for startAngles
if 1==1
    figure(23232);
    clf;
    hold on;
    grid on;

    % Fit the data?
    xdata = log(offsets_to_try);
    ydata = startAngle_intercepts;

    % Fit equation of form
    % y = x2*a + x*b + c --> Y = [X^3 X^2 X 1]*[a b c]'
    X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];
    % X_vec = [xdata.^2 xdata ones(length(xdata(:,1)),1)];
    % X_vec = [xdata ones(length(xdata(:,1)),1)];
    fit_params = (X_vec'*X_vec)\(X_vec'*ydata);
    fprintf(1,'\n\nINTERCEPT FITTING for startAngles: \n')

    fprintf(1,'Fit parameters for y = x^3 + x^2 + x + 1 for y=log(0.5-slope) and x=log(-log(offset)):\n')
    disp(fit_params);

    y_pred = X_vec*fit_params;
    plot(xdata,y_pred,'r.-','MarkerSize',20,'LineWidth',3);
    plot(xdata,ydata, 'b.','MarkerSize',10,'LineWidth',3);
    xlabel('log(offsets)');
    ylabel('startAngle intercepts');

    fprintf(1,'Comparison of actual and fit intercepts:\n');
    for ith_intercept = 1:length(startAngle_intercepts)
        offset = offsets_to_try(ith_intercept,1);
        fit_intercept = fcn_INTERNAL_estimateStartAngleInterceptFromOffset(offset);
        fprintf(1,'Offset: %.9f, startAngleIntercept Actual: \t %.9f  Fit: \t%.9f Ratio: %.4f\n',offset, startAngle_intercepts(ith_intercept,1),fit_intercept, startAngle_intercepts(ith_intercept,1)/fit_intercept);
    end
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

%% fcn_INTERNAL_estimateLengthSlopeFromOffset
function slope = fcn_INTERNAL_estimateLengthSlopeFromOffset(offset)
% This function estimates the slope from the offset

% Define the x-data used for fitting
xdata = log(-log(offset));

% xdata = log(-log(offsets_to_try));
% ydata = log(0.5-slopes);


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

fit_params = [-0.023713568015240
   0.205326944457505
  -0.596356120454393
   0.588196591187278];

y_pred = X_vec*fit_params;

slope = 0.5-exp(y_pred);
end % Ends fcn_INTERNAL_estimateLengthSlopeFromOffset

%% fcn_INTERNAL_estimateLengthInterceptFromOffset
function lengthIntercept = fcn_INTERNAL_estimateLengthInterceptFromOffset(offset)
% Define the x-data used for fitting
xdata = log(offset);

% xdata = log(offsets_to_try);
% ydata = intercepts;


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

% These were fit from earlier
fit_params = [
    -0.000648129505156
  -0.018419651621482
   0.315399921979639
   0.849677809445654
   ];

y_pred = X_vec*fit_params;

lengthIntercept = y_pred;
end % Ends fcn_INTERNAL_estimateLengthInterceptFromOffset


%% fcn_INTERNAL_estimateSpiralLength
function estimated_length = fcn_INTERNAL_estimateSpiralLength(circle2_radius, offset)

% x_data = log(1-circle2_radii);
% y_data = log(lengths);

x_data = log(1-circle2_radius);

% Get the slope and intercept
fit_slope = fcn_INTERNAL_estimateLengthSlopeFromOffset(offset);
fit_intercept = fcn_INTERNAL_estimateLengthInterceptFromOffset(offset);

estimated_length = exp(fit_slope*x_data + fit_intercept);
end % Ends fcn_INTERNAL_estimateSpiralLength


%% fcn_INTERNAL_estimateStartAngleSlopeFromOffset
function startAngleSlope = fcn_INTERNAL_estimateStartAngleSlopeFromOffset(offset)
% This function estimates the slope from the offset

    % Define the x-data used for fitting
xdata = log(-log(offset));

% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

fit_params = [ 
    0.027224755607716
  -0.232029761011433
   0.662948799520068
  -1.142906989473347];

y_pred = X_vec*fit_params;

startAngleSlope = y_pred;
end % Ends fcn_INTERNAL_estimateStartAngleSlopeFromOffset

%% fcn_INTERNAL_estimateStartAngleInterceptFromOffset
function startAngleIntercept = fcn_INTERNAL_estimateStartAngleInterceptFromOffset(offset)
% Define the x-data used for fitting
xdata = log(offset);

% xdata = log(offsets_to_try);
% ydata = startAngle_intercepts;


% Uses equation of form
% y = x2*a + x*b + c --> Y = [X^2 X 1]*[a b c]'
X_vec = [xdata.^3 xdata.^2 xdata ones(length(xdata(:,1)),1)];

% These were fit from earlier
fit_params = [
    -0.000805873132752
  -0.022273590488785
   0.285242206383266
   0.081271471353003
   ];

y_pred = X_vec*fit_params;

startAngleIntercept = y_pred;
end % Ends fcn_INTERNAL_estimateStartAngleInterceptFromOffset


%% fcn_INTERNAL_estimateSpiralStartAngle
function estimated_startAngle = fcn_INTERNAL_estimateSpiralStartAngle(circle2_radius, offset)

% %%%%%
% % Fit the start angle
% % Check the linear regression
% x_data = log(1-circle2_radii);
% y_data = log(start_angles);
%
% % Fit equation of form
% % y = x*m+b --> Y = [X 1]*[m b]
% X_vec = [x_data ones(length(x_data(:,1)),1)];
% fit_params = (X_vec'*X_vec)\(X_vec'*y_data);
% startAngle_slopes(ith_offset,1) = fit_params(1,1);
% startAngle_intercepts(ith_offset,1) = fit_params(2,1);
%
% y_pred = [x_data ones(length(x_data(:,1)),1)]*fit_params;
%

% x_data = log(1-circle2_radii);
% y_data = log(start_angles);

x_data = log(1-circle2_radius);

% Get the slope and intercept
fit_slope = fcn_INTERNAL_estimateStartAngleSlopeFromOffset(offset);
fit_intercept = fcn_INTERNAL_estimateStartAngleInterceptFromOffset(offset);

estimated_startAngle = exp(fit_slope*x_data + fit_intercept);
end % Ends fcn_INTERNAL_estimateSpiralStartAngle


