% Algorithm
%
% The algorithm takes the first two points (1st and 2nd point) and fits a
% straight line and an arc. Then, it calculates the RMSEs (root mean square
% error) for the straight line and arc fit and store the values in two
% different matrices (one is for RMSE of the straight line fit and the
% other one is for RMSE of the arc fit).
%
% Similarly, the algorithm takes the first three points (1st, 2nd and 3rd)
% and continues the process and store the RMSEs in their respective
% matrices. The algorithm takes the consecutive points and follows the same
% process for all the points. And store the RMSEs of all the combinations
% in their respective matrices
% 
% After determining the RMSE matrices for straight line and arc fit, the
% segment size is calculated.
%
% The optimal segment sizes from the current point to fit either an arc
% or a straight line are computed by adding the number of RMSEs that
% are less than error tolerance and greater than zero.
%
% Hence two different segment sizes are obtained, one by doing the above
% operation on RMSE matrix of the straight line fit and other by doing the
% same operation on RMSE matrix of the arc fit. Based on length of the
% segment sizes, the right fit is chosen.
%
% The greater value is chosen as the segment size and the fit that gave the
% greater value is chosen as the right fit. 
%
% The currentPoint is updated after fitting the each segment based on the
% segment size.
%
% This algorithm continues until the currentPoint is less than the length
% of the data.
%
% Example (Demo):
%
% Let's say there are 165 points of data and the error tolerance chosen by
% the user is 0.2. currentPoint = 1 and errorTolerance = 0.2
%
% Firstly, a for loop iterates from currentPoint = 1 to length of the data,
% i.e. 165.
%
% 1st iteration: 1st point --> Fits the point using an arc and a straight
% line. --> Calculates the RMSE and stores in the arcFit_RMSE and
% lineFit_RMSE matrices, respectively.
%
% 2nd iteration: 1st point and 2nd Point --> Fit the points using an arc
% and a straight line. --> Calculates the RMSE and stores in the
% arcFit_RMSE and lineFit_RMSE matrices, respectively.
%
% 3rd iteration: 1st, 2nd and 3rd Point --> Fit the points using an arc and
% a straight line. --> Calculates the RMSE and stores in the arcFit_RMSE
% and lineFit_RMSE matrices, respectively.
%
% 165th(last) iteration: 1st, 2nd, 3rd, 4th.....165th Point --> Fit the
% points using an arc and a straight line. --> Calculates the RMSE and
% stores in the arcFit_RMSE and lineFit_RMSE matrices, respectively.
%
% arcSegment_size and lineSegment_size are computed by summing the elements
% of arcFit_RMSE and lineFit_RMSE that are less than error tolerance and
% greater than zero.
%
% if arcSegment_size is greater than lineSegment_size, the arcSegment_size
% is chosen as the segment size and arc fit is used to fit that particular
% segment and viceversa.
%
% Let's assume arcSegment_size = 3 and lineSegment_size = 2, therefore, the
% segmentSize is chosen as 3 and arc fit is used to fit that segment.
%
% After determining the segment size, the currentPoint is updated to
% currentPoint = currentPoint + (segmentSize - 1)
%
% Therefore, the currentPoint = 3 (in this iteration).
%
% Similarly, now the loop iterates from 3 (currentPoint) to 165.
%
% 1st iteration: 3rd point --> Fits the point using an arc and a straight
% line. --> Calculates the RMSE and stores in the arcFit_RMSE and
% lineFit_RMSE matrices, respectively.
%
% 2nd iteration: 3rd and 4th Point --> Fit the points using an arc and a
% straight line. --> Calculates the RMSE and stores in the arcFit_RMSE and
% lineFit_RMSE matrices, respectively.
%
% 3rd iteration: 3rd, 4th and 5th Point --> Fit the points using an arc and
% a straight line. --> Calculates the RMSE and stores in the arcFit_RMSE
% and lineFit_RMSE matrices, respectively.
%
% 162th(last) iteration: 3rd, 4th, 5th.....165th Point --> Fit the points
% using an arc and a straight line. --> Calculates the RMSE and stores in
% the arcFit_RMSE and lineFit_RMSE matrices, respectively.
%
% arcSegment_size and lineSegment_size are computed by summing the elements
% of arcFit_RMSE and lineFit_RMSE that are less than error tolerance and
% greater than zero.
%
% if lineSegment_size is greater than arcSegment_size, the lineSegment_size
% is chosen as the segment size and arc fit is used to fit that particular
% segment and viceversa.
%
% Let's assume arcSegment_size = 3 and lineSegment_size = 5, therefore, the
% segmentSize is chosen as 5 and straight line fit is used to fit that segment.
%
% After determining the segment size, the currentPoint is updated to
% currentPoint = currentPoint + (segmentSize - 1)
%
% Therefore, the currentPoint = 7 (in this iteration).
%
% The algorithm stops when the currentPoint > 165.  



clc
close all

% Generating random data 
scenario = drivingScenario;
roadCenters = [0 0 0; 50 0 0; 50 50 0; 0 50 0; -50 50 0; -50 0 0; 0 0 0];
roadWidth = 3.5;
road(scenario, roadCenters, roadWidth);

% Extract only the x,y coordinates of the data
edge_coordinates_cells = roadBoundaries(scenario);
edge1 = edge_coordinates_cells{1};
edge2 = edge_coordinates_cells{2};

% Add a little noise to the data
edge1 = edge1 + 0.1 * randn(size(edge1));
edge2 = edge2 + 0.1 * randn(size(edge2));

% Plot the noisy data
f1 = figure('Color',[1,1,1]);
hold on;
p1 = plot(edge1(:, 1), edge1(:, 2), '.',DisplayName='Edge 1',MarkerSize=12);
% p2 = plot(edge2(:, 1), edge2(:, 2), '.',DisplayName='Edge 2');
title('Noisy Data Fitted using an Arc or a Straight Line');
xlabel('X');
ylabel('Y');

grid on
box on
legend

% Error tolerance for root mean square error. 
errorTolerance = 0.2; 

% lineFit_RMSE = [(1:size(edge1,1))', zeros(size(edge1,1),1)];
% arcFit_RMSE = [(1:size(edge1,1))', zeros(size(edge1,1),1)];

num_arcs = 0;
num_lines = 0;

fittedEdge = [];
currentPoint = 1;

% This loop runs until currentPoint is greater than the length of the edge1
while currentPoint < size(edge1,1)

    % The lineFit_RMSE and arcFit_RMSE are preallocated for the speed. 
    lineFit_RMSE = [(1:size(edge1,1))', zeros(size(edge1,1),1)];
    arcFit_RMSE = [(1:size(edge1,1))', zeros(size(edge1,1),1)];

    % This for loop calculates the RMS error from the current point to the
    % last point of the data set. 
    for i = currentPoint:size(edge1,1)

        noisyData = [edge1(currentPoint:i,1), edge1(currentPoint:i,2)];

        % Straight Line
        [fittedLine,fittedParametersLine] = fcn_geometry_fitStraightLine(noisyData);

        % RMS error when a straight line is used to fit the points
        RMSerrorStraightLine = sqrt(mean(sum((noisyData-fittedLine).^2,2))); % Calculating the RMS error for Line Fit

        lineFit_RMSE(i,2) = RMSerrorStraightLine;

        % Arc
        initial_center_xCoor_guess = mean(noisyData(:,1));
        initial_center_yCoor_guess = mean(noisyData(:,2));
        initial_radius_guess = 1;
        initial_startAngle_guess = 0;
        initial_endAngle_guess = pi/2;

        initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

        [fittedArc,fittedParametersArc] = fcn_geometry_fitArc(noisyData, initialGuess);
        
        % RMS error when an arc is used to fit the points
        RMSerrorArc = sqrt(mean(sum((noisyData-fittedArc).^2,2))); % Calculating the RMS error for an Arc Fit

        arcFit_RMSE(i,2) = RMSerrorArc;

    end

       % arcFit_RMSE = round(arcFit_RMSE,9);
       % lineFit_RMSE = round(lineFit_RMSE,9);
    
    % The optimal segment sizes from the current point to fit either an arc
    % or a straight line are computed by adding the number of RMSEs that
    % are less than error tolerance and greater than zero.
    arcSegment_size = sum(arcFit_RMSE(:,2) < errorTolerance & arcFit_RMSE(:,2) > 0); 
    lineSegment_size = sum(lineFit_RMSE(:,2) < errorTolerance & lineFit_RMSE(:,2) > 0); 

    % arcSegment_size = sum(arcFit_RMSE(:,2) < errorTolerance); 
    % lineSegment_size = sum(lineFit_RMSE(:,2) < errorTolerance);

    % Based on length of the segment sizes, the right fit is chosen. This
    % conditional statement fits an arc, if the arcSegment_size is bigger
    % than the lineSegment_size and viceversa.
    if arcSegment_size >= lineSegment_size
        
        segmentSize = arcSegment_size;
        % segmentSize
        endIndex = min(currentPoint+(segmentSize-1), size(edge1,1));
        noisyData = [edge1(currentPoint:endIndex,1), edge1(currentPoint:endIndex,2)];

        initial_center_xCoor_guess = mean(noisyData(:,1));
        initial_center_yCoor_guess = mean(noisyData(:,2));
        initial_radius_guess = 0.1;
        initial_startAngle_guess = 0;
        initial_endAngle_guess = pi/2;

        initialGuess = [initial_center_xCoor_guess,initial_center_yCoor_guess,initial_radius_guess,initial_startAngle_guess, initial_endAngle_guess];

        [fittedPoints,fittedParametersArc] = fcn_geometry_fitArc(noisyData, initialGuess, 1234);
        
        % Number of arcs used to fit the data are computed here. 
        num_arcs = num_arcs + 1;

    else

        segmentSize = lineSegment_size;
        % segmentSize
        endIndex = min(currentPoint+(segmentSize-1), size(edge1,1));
        noisyData = [edge1(currentPoint:endIndex,1), edge1(currentPoint:endIndex,2)];

        [fittedPoints,fittedParametersLine] = fcn_geometry_fitStraightLine(noisyData, 1234);
        
        % Number of lines used to fit the data are computed here. 
        num_lines = num_lines + 1;

    end
    
    % The fitted points are stored in this variable
    fittedEdge(currentPoint:endIndex,:) = fittedPoints;
    
    % The current point is updated after fitting the each segment 
    currentPoint = currentPoint+(segmentSize-1);

    % if 0 == currentPoint 
    %     currentPoint = 2;
    % end

end

figure(f1)
plot(fittedEdge(:,1),fittedEdge(:,2),'-r',LineWidth=1.5,DisplayName='Fitted Edge');






