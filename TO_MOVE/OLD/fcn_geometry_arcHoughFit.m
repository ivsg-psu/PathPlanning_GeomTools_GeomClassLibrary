function [fittedParameters, agreementIndices] = fcn_geometry_arcHoughFit(inputPoints, tolerance, varargin)
% fcn_geometry_arcHoughFit
%
% This function takes the input points and tolerance as the input and
% outputs the fitted parameters and agreement indices. 
% 
% [fittedParameters, agreementIndices] = fcn_geometry_circleHoughFit(inputPoints, tolerance)
% 
% INPUTS:
% 
% inputPoints: a  Nx2 vector where N is the number of points, but at
%              least 2.
%
% tolerance: The indices that are in the "tolerance" limit are stored.
%
% OUTPUTS:
%
% fittedParameters: The fitted parameters of the circle are stored here. 
% The size of the matrix is size(nchoosek(1:N,3),1)x2. 
%
% fittedParameters = [arcCenter, arcRadius, theta_start, theta_end]
%
% agreementIndices: The indices of the points that are within the tolerance
% limit are stored in this matrix. The size of this matrix is
% size(nchoosek(1:N,3), 1)xN. 
%
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%      
% See the script: script_test_fcn_geometry_arcHoughFit
% for a full test suite.
%
% This function was written on 2023_12_25 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu 

% Revision history:
% 2023_12_25 
% -- wrote the code: A. Batchu

flag_do_plots = 0; % Flag to plot
flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 34838;
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

if flag_check_inputs
    % Are there the right number of inputs?
    narginchk(3,4);

end

% Does user want to show the plots?
if 3 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1; 
    end
else
    if flag_do_debug
        flag_do_plots = 1;
    end
end

%% Main Code starts from here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total number of points
N_points = size(inputPoints,1);

% All possible 3-point combinations
combos_paired = nchoosek(1:N_points,3);
N_combos = size(combos_paired,1);

% Pre-allocation of fittedParameters and agreementIndices for saving
% computation time
fittedParameters = zeros(N_combos,5);
agreementIndices = zeros(N_combos, N_points);

for ith_combo = 1:N_combos
    
    % Inputs for fcn_geometry_arcCenterFrom3Points
    point1 = inputPoints(combos_paired(ith_combo,1),:);
    point2 = inputPoints(combos_paired(ith_combo,2),:);
    point3 = inputPoints(combos_paired(ith_combo,3),:);

    % fcn_geometry_arcCenterFrom3Points determines the fitted parameters
    [arcCenter, arcRadius, theta_start, theta_end] = fcn_geometry_arcCenterFrom3Points(point1, point2, point3, fig_num);
    % fitted parameters are stored in "fittedParameters" matrix
    fittedParameters(ith_combo,:) = [arcCenter, arcRadius, theta_start, theta_end];

    % Finding Agreement Indices

    % Distance of all the input points from the center
    distance_inputpoints_point1 = sum((inputPoints - point1).^2,2).^0.5;
    distance_inputpoints_point3 = sum((inputPoints - point3).^2,2).^0.5;

    % Root-mean-square error of distance_inputpoints_point1 and 
    % distance_inputpoints_point3 to find indices in agreement 
    RMS_error = sqrt(mean(sum((distance_inputpoints_point3-distance_inputpoints_point1).^2,2)));

    % Indices in agreement
    indicies_in_agreement = (RMS_error < tolerance)';
    
    % indices in agreement found in each iteration are stored in
    % "agreementIndices" matrix
    agreementIndices(ith_combo,:) = indicies_in_agreement;

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

end
if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

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