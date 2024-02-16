function [breakPointsCell, breakPointClosePairs] = fcn_geometry_findBreakpoints(domains)

% test script: script_test_fcn_geometry_findBreakpoints

% This function is still under testing. 
% The instructions will be updated after testing the function

N_houghDomains = size(domains,2) - 1;

% Empty breakPoints domain structure
breakPoints.firstBreakPoint = [nan nan];
breakPoints.lastBreakPoint = [nan nan];
breakPoints.fitType = 'empty';
breakPoints.fitParameters = nan;

% Create a cell array to save the structure of each fit
breakPointsCell = cell(1,N_houghDomains);

% This loop saves the first and last points of each fit domain as the first
% and last break points in breakPointsCell cell array. This cell array also
% stores the fit type and parameters

breakPointsMatrix = zeros(2*N_houghDomains, 2);

for i = 1:N_houghDomains

    breakPointsCell{i} = breakPoints;
    breakPointsCell{i}.firstBreakPoint = domains{i}.points_in_domain(1,:);
    breakPointsCell{i}.lastBreakPoint = domains{i}.points_in_domain(end,:);
    breakPointsCell{i}.fitType = domains{i}.best_fit_type;
    breakPointsCell{i}.fitParameters = domains{i}.best_fit_parameters;

    % Saving the break points in an array to find the closest break points
    breakPointsMatrix(2*(i-1)+1,:) = breakPointsCell{i}.firstBreakPoint;
    breakPointsMatrix(2*i,:) = breakPointsCell{i}.lastBreakPoint;

end

% All possible 2 point combinations 
combos = nchoosek(1:size(breakPointsMatrix,1), 2);

% Find distances between the break points
distance_btw_breakpoints = (sum((breakPointsMatrix(combos(:,2),:) - breakPointsMatrix(combos(:,1),:)).^2,2)).^0.5;

% If the distance between two break points are less than or equal to this
% tolerance, they are considered as the break point close pairs
tolerance = 0.5;

% Finding the indices of break point close pairs
idx_breakpoint_close_pairs = combos(distance_btw_breakpoints <= tolerance,:);

% Reshaping the indices as a column matrix
idx_breakpoint_close_pairs_column = reshape(idx_breakpoint_close_pairs.', [], 1);

% Finding the break point close pairs 
breakPointClosePairs = breakPointsMatrix(idx_breakpoint_close_pairs_column,:);

end


