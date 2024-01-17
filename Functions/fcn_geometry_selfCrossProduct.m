function [cross_products, err] ...
    = ...
    fcn_geometry_selfCrossProduct(...
    path,...
    varargin)
% fcn_geometry_selfCrossProduct
% finds self cross product of a path, e.g. whether at each internal point
% of a path, the path bends next to the left (positive cross product) or to
% the right (negative cross product).
%
% FORMAT:
%
% [cross_products] ...
%     = ...
%     fcn_geometry_selfCrossProduct(...
%     path,...
%     (fig_num))
%
% INPUTS:
%
%      path: an [N x 2] vector of X,Y data for the path, with N>=3.
%
%      (OPTIONAL INPUTS)
%  
%      fig_num: a figure number to plot results.
%
% OUTPUTS:
%
%      cross_products: a [(N-2) x 1] vector of cross products calculated by
%      crossing segment 1 to segment 2, segment 2 to segment 3. If there N
%      points in the original path, there will be N-1 segments and
%      therefore N-2 cross product results. The cross products are
%      organized such that the first cross-product is between segments 1
%      and 2, the second is between 2 and 3, etc.
% 
%      err: returns 0 (good) if all the cross-products are well defined, 1
%      if there was an error, e.g. any corners are ill-defined such as
%      stright lines, zero lengths, or segments bent back on themselves.
%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_geometry_findTangentPointsFromPointToCircle
%
% EXAMPLES:
%
% See the script: script_test_fcn_geometry_selfCrossProduct
% for a full test suite.
%
% This function was written on 2021_04_25 by S. Brennan
% by modifying: fcn_geometry_findTangentPointsTwoCircles
% Questions or comments? sbrennan@psu.edu 
%
% See: http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm for
% the mathematics being implemented here.

% Revision History:
% 2021-04-25
% -- First write of the code
% 2024_01_08 - S. Brennan
% -- fixed bug with cross function call to force it to cross column-wise

%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments.')
    end
    
    % Check the path input
    fcn_DebugTools_checkInputsToFunctions(...
        path, '2column_of_numbers',[3 4]);
    
end

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_do_plot = 1;
else
    if flag_do_debug
        fig = figure;
        fig_num = fig.Number;
        flag_do_plot = 1;
    end
end


%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncross = length(path(:,1))-2;  % The number of cross products to calculate
zero_col = zeros(Ncross,1);

%% Initialize vectors
first_vectors = diff(path(1:end-1,:));
second_vectors = diff(path(2:end,:));

% Calculate the cross product
cross_result = cross([first_vectors zero_col],[second_vectors zero_col],2);
cross_products = cross_result(:,3);

% Check for errors?
err = 0;
if nargout > 1    
    if any(cross_products == 0)
        err = 1;
    end
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
if flag_do_plot
    figure(fig_num);
    hold on;
    axis equal;
    grid on; grid minor;
    
    plot(path(:,1),path(:,2),'.-');
    
    for  i=1:Ncross
        text(path(i+1,1),path(i+1,2),sprintf('Cross is: %.2f',cross_products(i)));
    end
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end

end % Ends function


