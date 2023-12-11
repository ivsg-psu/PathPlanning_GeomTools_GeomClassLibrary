% script_test_fcn_geometry_fillLineTestPoints
% Exercises the function: fcn_geometry_fillLineTestPoints
% Revision history:
% 2023_12_05
% -- wrote the code

close all;
clc;


%% Fill test data - 3 segments
fig_num = 23;
seed_points = [2 3; 4 5; 7 0; 9 5];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,fig_num);


%% Calculate vectors from from all points to each other 
% This is a NxN calculation
input_points = test_points;

N_points = length(input_points(:,1));
points_repeated = repelem(input_points,N_points,1);
points_stacked  = repmat(input_points,N_points,1);

figure(38383);
clf;
hold on;
grid on;
axis equal;

for ith_point = 1:N_points-1
    for jth_point = ith_point+1:N_points
        projection_vector = input_points(jth_point,:) - input_points(ith_point,:);
        unit_projection = fcn_INTERNAL_calcUnitVector(projection_vector);
        plot(unit_projection(:,1),unit_projection(:,2),'b.','MarkerSize',20);
    end
end

%% Solving using permutations
combos_paired = nchoosek(1:N_points,2);
projection_vectors = input_points(combos_paired(:,2),:) - input_points(combos_paired(:,1),:);
unit_projection_vectors = fcn_INTERNAL_calcUnitVector(projection_vectors);
angles = atan2(unit_projection_vectors(:,2),unit_projection_vectors(:,1));


figure(38383);
clf
hold on;

hist(angles,360)

% Create weighted histogram
agreements = zeros(size(angles));
threshold = 0.06; % Units are meters

for ith_vector = 1:length(projection_vectors(:,1))
    unit_orthogonal_vector = unit_projection_vectors(ith_vector,:)*[0 1; -1 0];
    base_projection_vectors = input_points - input_points(combos_paired(ith_vector,1),:);
    lateral_distances = sum(unit_orthogonal_vector.*base_projection_vectors,2);

    agreements(ith_vector) = sum(abs(lateral_distances)<threshold);

end

figure(2334);
clf
hold on;
grid on;
plot(angles,agreements,'k.');


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

function unit_projection = fcn_INTERNAL_calcUnitVector(projection_vector)
vector_length = sum(projection_vector.^2,2).^0.5;
unit_projection = projection_vector./vector_length;
end
