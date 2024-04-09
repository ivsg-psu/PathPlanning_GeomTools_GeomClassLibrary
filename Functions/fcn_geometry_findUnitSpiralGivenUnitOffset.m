function [spiralLength, spiralEndAngleInRadians] = fcn_geometry_findUnitSpiralGivenUnitOffset(varargin)
%% fcn_geometry_findUnitSpiralGivenUnitOffset
% Given a line segment that is along the -x axis increasing from negative
% to positive, and ending at the origin, and given a circle of unit radius
% that this line will "spiral" into in the positive direction, and given
% the y-distance offset from the bottom of this circle above the x-axis,
% finds the spiral that would join the line segment to the unit circle
%
% The method to solve this is to use repeated spiral generations to
% approximate the relationship between the y-offset, the spiral length, and
% the spiral ending orientation. This is then used as a look-up fit to find
% the relationship from y-offset as an input to the other terms. The
% spirals are generated from fcn_geometry_extractXYfromSTSpiral.
%
% FORMAT: 
%
%       [x_arc,y_arc] = fcn_geometry_findUnitSpiralGivenUnitOffset(offset, varargin)
%
% INPUTS:
%
%      (OPTIONAL INPUTS): 
%
%      offset : the y-axis difference between the incoming line segment and
%      the bottom of the circle used as a target to design the spiral. 
% 
%      NOTE: If left empty, causes the the unit spiral calculations to be
%      performed. These calculations are persistent to the function and as
%      well are saved to file to avoid calculation in the future, as this
%      process is very slow.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%       spiralLength: the length of the spiral from the end of the line
%       segment to the start of the arc
%
%       spiralEndAngleInRadians: the ending angle of the spiral. NOTE: this
%       value is always equal to the spiralLength/2.
%
% DEPENDENCIES:
%
%      None
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_geometry_findUnitSpiralGivenUnitOffset.m for
%       a full test suite.
%
% This function was written by S. Brennan on 2024_04_06
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_04_06 - S. Brennan
% -- wrote the code

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==2 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_PARSEXODR_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_PARSEXODR_FLAG_CHECK_INPUTS");
    MATLABFLAG_PARSEXODR_FLAG_DO_DEBUG = getenv("MATLABFLAG_PARSEXODR_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_PARSEXODR_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_PARSEXODR_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_PARSEXODR_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_PARSEXODR_FLAG_CHECK_INPUTS);
    end
end

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(0,2);
    end
end


% Does user want to specify flag_initialize_values?
flag_initialize_values = 0;
if (0 ==nargin)
    flag_initialize_values = 1;
    offset = 0;
elseif (1<=nargin)
    temp = varargin{1};
    if isempty(temp)
        flag_initialize_values = 1;
    else
        offset = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (2<= nargin)
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end

%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if flag_initialize_values
    disp('Preinitialization is not yet coded.');
else
    spiralLength = fcn_INTERNAL_findLengthFromOffset(offset);
    spiralEndAngleInRadians = spiralLength;

end

%% THE FOLLOWING CODE IS TO PRECALCULATE THE RESULTS FOR SPEED
% It does not yet work, but is 99% done

% persistent unit_spiral_offsets unit_spiral_end_angles_radians unit_spiral_lengths unit_spiral_end_XYs
% 
% % Fill in initial parameters - only does this once per MATLAB instance
% if 1==flag_initialize_values || isempty(unit_spiral_offsets) || isempty(unit_spiral_end_angles_radians) || isempty(unit_spiral_lengths) || isempty(unit_spiral_end_XYs)
% 
%     % Check to see if data was loaded earlier
%     mat_filename = fullfile(cd,'Data','UnitSpiralFittingData.mat'); 
% 
%     if 1==flag_initialize_values
%         flag_load_all_data = 1;
%     else
% 
%         flag_load_all_data = 0; % Default value
% 
%         % Does the file exist?
%         if exist(mat_filename,'file')
%             load(mat_filename,'unit_spiral_offsets','unit_spiral_end_angles_radians','unit_spiral_lengths','unit_spiral_end_XYs');
%         else
%             % File does not exist - need to load it
%             flag_load_all_data = 1;
%         end
%     end
% 
% 
% 
%     if flag_load_all_data
% 
%         % Load the data from scratch
%         [unit_spiral_offsets,unit_spiral_end_angles_radians, unit_spiral_lengths, unit_spiral_end_XYs] = fcn_INTERNAL_calculateUnitSpiralData(fig_num);
% 
%         % % Grab the file's date of creation
%         % st = dbstack;
%         % this_function = st(1).file;
%         % file_info = dir(which(this_function));
%         % permanent_file_date = file_info.date;
% 
%         % Save the results
%         save(mat_filename,'unit_spiral_offsets','unit_spiral_end_angles_radians','unit_spiral_lengths','unit_spiral_end_XYs');
% 
%     end
% 
% end



%% Plot the results (for debugging)?
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
if flag_do_plots
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    % % Determine the hold state of the figure
    % holdState = ishold;
    % % Turn on hold for the purposes of adding the clothoid segment
    % hold on
    % 
    % grid on;
    % axis equal
    % 
    % 
    % % Grab colors to use
    % try
    %     color_ordering = orderedcolors('gem12');
    % catch
    %     color_ordering = colororder;
    % end
    % N_colors = length(color_ordering(:,1)); 
    % 
    % % Set up labels and title
    % xlabel('X (m)')
    % ylabel('Y (m)')
    % 
    % % Plot the results as dots
    % for ith_column = 1:length(x_spiral(1,:))
    %     % Get current color
    %     current_color = color_ordering(mod(ith_column,N_colors)+1,:); %#ok<NASGU>
    % 
    %     % Plot results as dots
    %     plot(x_spiral(:,ith_column),y_spiral(:,ith_column),'.-.','MarkerSize',30); % ,'Color',current_color);
    % end
    % 
    % 
    % % Restore the original hold state of the plot as necessary
    % if 0 == holdState
    %     hold off
    % end

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end


end % Ends check if plotting

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

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

%% fcn_INTERNAL_findLengthFromOffset
function sprialLength = fcn_INTERNAL_findLengthFromOffset(offset)
function_to_optimize = @(x)fcn_INTERNAL_calcUnitSpiralOffsetError(x,offset);
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% sprialLength = fminsearch(function_to_optimize,1,options);
sprialLength = fminsearch(function_to_optimize,1);

end % Ends fcn_INTERNAL_findLengthFromOffset

%% fcn_INTERNAL_findLengthFromOffset
function error = fcn_INTERNAL_calcUnitSpiralOffsetError(length_spiral,offset)

% PERFORM UNIT CIRCLE ANALYSIS
% Initial heading, position, and curvature of the spiral
h0 = 0;
x0 = 0;
y0 = 0;
K0 = 0;
Kf = 1; % Final curvature forced to have a radius of 1

% Call the function
[x_arc,y_arc] = fcn_geometry_extractXYfromSTSpiral(length_spiral,length_spiral,h0,x0,y0,K0,Kf,-1);

unit_spiral_end_angles_radians = length_spiral/2;
unit_tangent = [cos(unit_spiral_end_angles_radians) sin(unit_spiral_end_angles_radians)];

unit_orthogonal = unit_tangent*[0 1; -1 0];
circle_center = (1/Kf)*unit_orthogonal + [x_arc(end) y_arc(end)];

actual_offset            = circle_center(1,2) - (1/Kf);
error = abs(offset-actual_offset);


end % Ends fcn_INTERNAL_findLengthFromOffset

%% fcn_INTERNAL_calculateUnitSpiralData
function [unit_spiral_offsets,unit_spiral_end_angles_radians, unit_spiral_lengths, unit_spiral_end_XYs] = fcn_INTERNAL_calculateUnitSpiralData(fig_num) %#ok<DEFNU>

% PREP THE FIGURE?
% Plot the circle
if ~isempty(fig_num)
    figure(fig_num);
    clf;
    subplot(2,2,1); axis equal; grid on; hold on; xlabel('Unit x [m]'),ylabel('Unit y [m]');
    subplot(2,2,2); axis equal; grid on; hold on; xlabel('Unit x [m]'),ylabel('Unit y [m]');
    subplot(2,2,3); axis equal; grid on; hold on; xlabel('Unit x [m]'),ylabel('Unit y [m]');
    subplot(2,2,4); axis equal; grid on; hold on; xlabel('Unit x [m]'),ylabel('Unit y [m]');
    
end

% PERFORM UNIT CIRCLE ANALYSIS
% Initial heading, position, and curvature of the spiral
h0 = 0;
x0 = 0;
y0 = 0;
K0 = 0;
Kf = 1; % Final curvature forced to have a radius of 1

N_points = 40;
unit_spiral_offsets            = zeros(N_points,1);
unit_spiral_end_angles_radians = zeros(N_points,1);
unit_spiral_lengths            = linspace(0.01,6,N_points)';
unit_spiral_end_XYs            = zeros(N_points,2);


for length_index = 1:N_points
    if 0==mod(length_index,5)
        fprintf(1,'Testing length %.0d of %.0d\n',length_index,N_points);
    end
    L0 = unit_spiral_lengths(length_index);

    if length_index == N_points
        s  = 0:0.1:L0;
    else
        s  = L0;
    end

    % Call the function
    [x_arc,y_arc] = fcn_geometry_extractXYfromSTSpiral(s,L0,h0,x0,y0,K0,Kf,-1);

    unit_spiral_end_XYs(length_index,:) = [x_arc(end) y_arc(end)];

    % NOTE: final angle is exactly half the length for euler spiral. To
    % prove this, uncomment the if statement below
    if 1==0
        % NUMERICALLY find the center of the circle tangent at the end of the
        % spiral Find the unit vector (need to do this analytically!)
        s_tangent = [0.99999999 1]'*L0;
        [x_tangent,y_tangent] = fcn_geometry_extractXYfromSTSpiral(s_tangent,L0,h0,x0,y0,K0,Kf);
        unit_tangent = fcn_geometry_calcUnitVector([diff(x_tangent) diff(y_tangent)]);
        unit_spiral_end_angles_radians(length_index,:) = atan2(unit_tangent(2),unit_tangent(1));
    else
        unit_spiral_end_angles_radians(length_index,:) = L0/2;
        unit_tangent = [cos(unit_spiral_end_angles_radians(length_index,:)) sin(unit_spiral_end_angles_radians(length_index,:))];
    end

    unit_orthogonal = unit_tangent*[0 1; -1 0];
    circle_center = (1/Kf)*unit_orthogonal + [x_arc(end) y_arc(end)];

    unit_spiral_offsets(length_index,:)            = circle_center(1,2) - (1/Kf);


    % Plot the circle
    if 0==mod(length_index,100) && ~isempty(fig_num)
        figure(fig_num);
        subplot(2,2,1);
        plot(circle_center(:,1),circle_center(:,2),'r+');
        angles = (0:1:360)'*pi/180;
        XY_circle = (1/Kf)*[cos(angles) sin(angles)] + circle_center;
        plot(XY_circle(:,1),XY_circle(:,2),'r-');

        pause(0.01);
    end
end

figure(fig_num);
subplot(2,2,1);
plot([-1 0],[0 0],'k-');
plot(x_arc,y_arc,'b-');
plot(circle_center(:,1),circle_center(:,2),'r+');
angles = (0:1:360)'*pi/180;
XY_circle = (1/Kf)*[cos(angles) sin(angles)] + circle_center;
plot(XY_circle(:,1),XY_circle(:,2),'r-');

figure(fig_num);
subplot(2,2,2);
plot(unit_spiral_lengths, unit_spiral_offsets,'k.-');
xlabel('Spiral lengths [m]');
ylabel('Sprial offsets [m]');


figure(fig_num);
subplot(2,2,3);
plot(unit_spiral_lengths,unit_spiral_end_angles_radians,'k.-');
xlabel('Sprial lengths [m]');
ylabel('Spiral end angle [rad]');



% Fit the data
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'rat55' );  % Has about 5 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'poly9' );  % Has about 7 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'exp1' );  % Has about 50 cm error, very smooth
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'exp2' );  % Has about 12 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'fourier4' );  % Has about 8 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'gauss8' );  % Has about 12 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'gauss1' );  % Has about 50 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'power2' );  % Has about 12 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'power1' );  % Has about 15 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'sin8' );  % Has about 6 cm error
% [L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'spline' );  % Has about 2 cm error
[L0_curve, goodness] = fit(unit_spiral_offsets, unit_spiral_lengths, 'pchip' );  % Has about 2 cm error - almost zero through entire test

disp(goodness);

% Check the fit by running all the mid-points
unit_spiral_offsets_test = zeros(N_points,1);
unit_spiral_test_lengths = [0.001; unit_spiral_lengths(1:end-1,1)+diff(unit_spiral_lengths)/2];

% Use lengths to predict offsets
for length_index = 1:N_points
    if 0==mod(length_index,5)
        fprintf(1,'Testing length %.0d of %.0d\n',length_index,N_points);
    end
    L0 = unit_spiral_test_lengths(length_index);

    if length_index == N_points
        s  = 0:0.1:L0;
    else
        s  = L0;
    end

    % Call the function
    [x_arc,y_arc] = fcn_geometry_extractXYfromSTSpiral(s,L0,h0,x0,y0,K0,Kf,-1);

    unit_spiral_end_angles_radians_test = L0/2;
    unit_tangent = [cos(unit_spiral_end_angles_radians_test) sin(unit_spiral_end_angles_radians_test)];

    unit_orthogonal = unit_tangent*[0 1; -1 0];
    circle_center = (1/Kf)*unit_orthogonal + [x_arc(end) y_arc(end)];

    unit_spiral_offsets_test(length_index,:)            = circle_center(1,2) - (1/Kf);

end

% Check the fit using the offsets to now predict the lengths
unit_spiral_fit_lengths = L0_curve(unit_spiral_offsets_test);
length_errors = unit_spiral_test_lengths - unit_spiral_fit_lengths;

figure(fig_num);
subplot(2,2,4);
plot(unit_spiral_test_lengths,length_errors,'k.');
axis normal;
xlabel('Sprial length [m]');
ylabel('Error in spiral length prediction [m]');

fprintf(1,'Done calculating calibration data.\n');

% % Find the maximum offset we will be testing
% max_L0 = 6;
% max_s  = max_L0;
% [x_arc,y_arc] = fcn_geometry_extractXYfromSTSpiral(max_s,max_L0,h0,x0,y0,K0,Kf,-1);
% max_unit_spiral_end_angles_radians = max_L0/2;
% unit_tangent = [cos(max_unit_spiral_end_angles_radians) sin(max_unit_spiral_end_angles_radians)];
% unit_orthogonal = unit_tangent*[0 1; -1 0];
% circle_center = (1/Kf)*unit_orthogonal + [x_arc(end) y_arc(end)];
% max_unit_spiral_offset = circle_center(1,2) - (1/Kf);
% 
% % Check the approximations from offsets to lengths, using a LOT of points
% % from 0 all the way to maximum offset
% N_confirmation_points      = 10*N_points;
% true_unit_spiral_offsets   = linspace(0,max_unit_spiral_offset,N_confirmation_points)';
% L0_predicted               = zeros(N_confirmation_points,1);
% errors_in_length           = zeros(N_confirmation_points,1);
% 
% % Loop through each offset, using it to predict the length, then using
% % length to calculate actual offset, and then using actual offset to again
% % estimate length. Use the first and second length calculations to estimate
% % the fitting error.
% for offset_index = 1:N_confirmation_points
%     if 0==mod(offset_index,5)
%         fprintf(1,'Testing offset %.0d of %.0d\n',offset_index,N_confirmation_points);
%     end
%     offset_input = true_unit_spiral_offsets(offset_index);
% 
%     % Call the fitting function
%     L0_predicted(offset_index,1) = interp1(unit_spiral_offsets,unit_spiral_lengths,offset_input,'pchip');
% 
%     % Check the true answer
%     s  = L0;
%     [x_arc, y_arc] = fcn_geometry_extractXYfromSTSpiral(s,L0_predicted(offset_index,1),h0,x0,y0,K0,Kf,-1);
% 
%     % Calculate the actual offset
%     unit_spiral_end_angles_radians_test= L0_predicted(offset_index,1)/2;
%     unit_tangent = [cos(unit_spiral_end_angles_radians_test) sin(unit_spiral_end_angles_radians_test)];
%     unit_orthogonal = unit_tangent*[0 1; -1 0];
%     circle_center = (1/Kf)*unit_orthogonal + [x_arc(end) y_arc(end)];
%     true_offset  = circle_center(1,2) - (1/Kf);
% 
%     % Convert true offset into approximate true length
%     L0_approximate_true = interp1(unit_spiral_offsets,unit_spiral_lengths,true_offset,'pchip');
% 
%     % Find error
%     errors_in_length(offset_index,1) = L0_approximate_true - L0_predicted(offset_index,1);
% end
% 
% figure(fig_num);
% subplot(2,2,4);
% plot(true_unit_spiral_offsets,errors_in_length,'k.-');
% xlabel('Sprial offset input [m]');
% ylabel('Error in spiral length prediction [m]');
% axis normal;


end % Ends fcn_INTERNAL_calculateUnitSpiralData