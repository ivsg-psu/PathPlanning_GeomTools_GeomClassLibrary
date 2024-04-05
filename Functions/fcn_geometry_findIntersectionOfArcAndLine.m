function intersectionPoints = fcn_geometry_findIntersectionOfArcAndLine(endPointsCell, tole, fig_num)

lineSegment = 1;
regArc = 2; 

figure(fig_num)

% Location of S (Start Point of Line Segment) 
pointS = endPointsCell{lineSegment}.firstEndPoint;

% Location of E (End Point of Line Segement)
pointE = endPointsCell{lineSegment}.lastEndPoint;

% Location of Arc center
pointC = endPointsCell{regArc}.fitParameters(1,1:2);

% Plot Arc Center
plot(pointC(1,1), pointC(1,2), 'c.','MarkerSize',30)

% Calculate vector EC
vectorEC = pointC - pointE;

% Calculate vector ES
vectorES = pointS - pointE;

% vector_magnitude = sum(input_vector.^2,2).^0.5;
% unitVector_input_vector= input_vector./vector_magnitude;

% Plot the vectorEC in GREEN
quiver(pointE(1,1), pointE(1,2), vectorEC(1,1), vectorEC(1,2), 'Color', 'g', 'LineWidth', 4);

% Plot the vectorES in BLUE
quiver(pointE(1,1), pointE(1,2), vectorES(1,1), vectorES(1,2), 'Color', 'b', 'LineWidth', 4);

% Calculate the unit vector of vectorES
vectorES_magnitude = sum(vectorES.^2,2).^0.5;
unit_vectorES = vectorES./vectorES_magnitude;

% Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
quiver(pointE(1,1), pointE(1,2), unit_vectorES(1,1), unit_vectorES(1,2), 'Color', 'y', 'LineWidth', 3);

% Calculate unit orthogonal vector of vectorES
unit_orthogonal_vectorES = unit_vectorES*[0 1; -1 0];

% Plot the unit orthogonal vector ES (unit_orthogonal_vectorES) in RED
quiver(pointE(1,1), pointE(1,2), unit_orthogonal_vectorES(1,1), unit_orthogonal_vectorES(1,2), 'Color', 'r', 'LineWidth', 3);

% Find the distance between the arc center and start point of the line
% segment by calculating the dot product of VectorEC and unit_orthogonal_vectorES
dist_btw_pointC_and_pointS = dot(vectorEC, unit_orthogonal_vectorES);

% Radius of the regression Arc
radiusArc = endPointsCell{regArc}.fitParameters(1,3);

if dist_btw_pointC_and_pointS > radiusArc
    disp('No Intersection')
end

% tole = 0.1;
% One intersection point, if the distance between arc center and start
% point of segment is equal to radius
if (abs(dist_btw_pointC_and_pointS) <= (radiusArc + tole)) && (abs(dist_btw_pointC_and_pointS) >= (radiusArc - tole))
    intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    disp(intersectionPoints)
    plot(intersectionPoints(1,1), intersectionPoints(1,2), '.', 'Color', 'r', 'MarkerSize',40)
end


% Two intersection points, if the distance between arc center and start
% point of segment is less than others
if abs(dist_btw_pointC_and_pointS) < (radiusArc - tole)

    % Center point of the two intersection points
    % centerPoint_of_intersectionPoints = vectorES + dot(vectorEC, unit_vectorES)*unit_vectorES;
    centerPoint_of_intersectionPoints = pointE + dot(vectorEC, unit_vectorES)*unit_vectorES;
    plot(centerPoint_of_intersectionPoints(1,1), centerPoint_of_intersectionPoints(1,2), '.', 'Color', 'g', 'MarkerSize',30);

    % Opposite side of right angle triangle: dist_btw_pointC_and_pointS
    dist_btw_centerPoint_and_pointC = dot(vectorEC, unit_orthogonal_vectorES);

    % Adjacent side of the right angle triangle
    dist_btw_centerPoint_and_intersectionPoints = abs((radiusArc^2 - dist_btw_centerPoint_and_pointC^2))^0.5;
    
    %
    intersectionPoint1 = centerPoint_of_intersectionPoints + dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoint2 = centerPoint_of_intersectionPoints - dist_btw_centerPoint_and_intersectionPoints*unit_vectorES; 

    intersectionPoints = [intersectionPoint1; intersectionPoint2]; 
    disp(intersectionPoints)
    plot(intersectionPoints(:,1), intersectionPoints(:,2), '.', 'Color', 'r', 'MarkerSize',40)
end
end