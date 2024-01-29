function unit_orthogonal_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, varargin)
% fcn_geometry_calcOrthogonalVector
% Finds the unit vectors orthogonal to an input list of input vectors. For
% 2D inputs, the orthogonal vectors are always positive direction. For N-D
% vectors, the orthogonal vectors are always orthogonal, but have random
% rotations around the axis of the input vectors.
%
% Format: 
% unit_orthogonal_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, (fig_num))
%
% INPUTS:
%      input_vectors: a list of Nxm vector where N is the number of vectors
%      that should be converted into unit-length orthogonal vector, and m
%      is the dimension of the vector (typically 2 or 3).
%
%      (OPTIONAL INPUTS)
% 
%      seed_points: a set of points, of Nxm, representing seed points for
%      calculating the orthogonal directions if the dimension is 3 or more.
%      Default is a random number.
% 
%      fig_num: a figure number to plot the results.
%
% OUTPUTS:
%
%      unit_orthogonal_vectors: the unit-length vectors orthogonal to each input vector.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_calcUnitVector
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_calcOrthogonalVector
% for a full test suite.
%
% This function was written on 2024_01_24 by S. Brennan
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2024_01_24 - S. Brennan
% -- wrote the code
% 2024_01_28 - S. Brennan
% -- added seed point inputs

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS");
    MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG = getenv("MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_GEOMETRY_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_GEOMETRY_FLAG_CHECK_INPUTS);
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

if 0==flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(1,3);

        % Check the projection_vector input to be length greater than or equal to 1
        fcn_DebugTools_checkInputsToFunctions(...
            input_vectors, '2or3column_of_numbers');

    end
end

% Does user want to specify fig_num?
seed_points = []; % Default is to have no seed pionts
if (2<= nargin) 
    temp = varargin{1};
    if ~isempty(temp)
        seed_points = temp;
    end
end

% Does user want to specify fig_num?
fig_num = []; % Default is to have no figure
flag_do_plots = 0;
if (0==flag_max_speed) && (3<= nargin) 
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        flag_do_plots = 1;
    end
end


%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many points do we have?
N_vectors = length(input_vectors(:,1));

% Calculate unit vectors
unit_input_vectors = fcn_geometry_calcUnitVector(input_vectors,-1);

% Find unit vectors orthogonal to point-to-point vectors
if length(input_vectors(1,:))==2 % 2D vectors
    unit_orthogonal_vectors = unit_input_vectors*[0 1; -1 0];

elseif length(input_vectors(1,:))==3 % 3D vectors

    if isempty(seed_points)
        % Generate random direction vectors
        random_directions = 2*rand(N_vectors,3) - ones(N_vectors,3);
    else
        random_directions = seed_points;
    end
    
    % Calculate the portion of the random direction aligned with the unit
    % vector. This is just the dot product
    distance_along_unit_vector = sum(random_directions.*unit_input_vectors,2);

    % Subtract off the aligned part. What is left is only the portion of
    % random directions orthogonal to the unit vectors
    random_orthogonal_vectors = random_directions - distance_along_unit_vector.*unit_input_vectors;

    % Convert to unit vectors
    unit_orthogonal_vectors = fcn_geometry_calcUnitVector(random_orthogonal_vectors,-1);

else
    error('unknown dimension encountered when calculating unit vectors');
end

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

    hold on;
    grid on;
    axis equal
    xlabel('X [m]')
    ylabel('Y [m]')

    % Plot the input vectors alongside the unit vectors
    if length(input_vectors(1,:))==2 % 2D vectors
        for ith_vector = 1:N_vectors
            h_plot = quiver(0,0,input_vectors(ith_vector,1),input_vectors(ith_vector,2),0,'-','LineWidth',3);
            plot_color = get(h_plot,'Color');
            quiver(0,0,unit_orthogonal_vectors(ith_vector,1),unit_orthogonal_vectors(ith_vector,2),0,'-','LineWidth',1,'Color',(plot_color+[1 1 1])/2,'MaxHeadSize',1);
        end

    else % 3D vectors
        zlabel('Z [m]')
        view(3);
        for ith_vector = 1:N_vectors
            h_plot = quiver3(0,0,0, input_vectors(ith_vector,1),input_vectors(ith_vector,2),input_vectors(ith_vector,3),0,'-','LineWidth',3);
            plot_color = get(h_plot,'Color');

            % Plot the results
            quiver3(input_vectors(ith_vector,1),input_vectors(ith_vector,2),input_vectors(ith_vector,3),...
                unit_orthogonal_vectors(ith_vector,1),unit_orthogonal_vectors(ith_vector,2),unit_orthogonal_vectors(ith_vector,3),0,'-','LineWidth',1,'Color',(plot_color+[1 1 1])/2,'MaxHeadSize',1);
        end
    end

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




