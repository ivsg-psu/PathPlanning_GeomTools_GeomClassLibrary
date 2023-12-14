% script_test_fcn_geometry_fillLineTestPoints
% Exercises the function: fcn_geometry_fillLineTestPoints

% Revision history:
% 2023_12_05
% -- wrote the code
% 2023_12_12
% -- updated for weighted Hough transform test

close all;
clc;


%% Fill test data - 3 segments
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 7 0; 9 5; 1 1; 13 14];
M = 10;
sigma = 0.2;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,fig_num);

%% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

%% Calculate vectors from from all points to each other 

lateral_match_threshold = 0.06; % Units are meters
threshold_max_points = 10;
input_points = test_points;
fig_points_plot = 12828;

% Main part of function starts here (HoughSegmentation)
N_points = length(input_points(:,1));


% Find all possible 2-point permutations
combos_paired = nchoosek(1:N_points,2);
N_combos = length(combos_paired(:,1));

projection_vectors = input_points(combos_paired(:,2),:) - input_points(combos_paired(:,1),:);
unit_projection_vectors = fcn_INTERNAL_calcUnitVector(projection_vectors);
angles = atan2(unit_projection_vectors(:,2),unit_projection_vectors(:,1));



% Plot the input points
figure(fig_points_plot);
clf;
hold on;
grid on;
axis equal;
plot(input_points(:,1),input_points(:,2),'k.','MarkerSize',20);



% Find the agreements between points and every single line fit
[agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(input_points,unit_projection_vectors, combos_paired, lateral_match_threshold);
agreements = sum(agreement_indicies,2);


% Find domains based on removing peaks one at a time

% Initialize the search process for domains by recording the counts of
% total agreement (agreements) and the points associated with the agreement
% counts

remaining_agreements = agreements;
remaining_points = input_points;
max_original_agreement = max(agreements);


domain_count = 0;
domains{1} = struct;

while(max(remaining_agreements)>threshold_max_points)
    domain_count = domain_count+1;

    % Plot the original agreements
    figure(fig_points_plot+1);
    if 1==domain_count
        clf
        hold on;
        grid on;
        tiledlayout('flow')

        % Create histogram of raw angles
        nexttile;
        histogram(angles,360)

    end    
    nexttile;
    plot(angles,remaining_agreements,'k.','MarkerSize',20);
    hold on;
    ylim([0 max_original_agreement*1.1]);


    % Pull out the maximum agreement
    [max_agreement,index_max_agreement] = max(remaining_agreements);
    indicies_in_agreement = agreement_indicies(index_max_agreement,:);
    indicies_of_points = find(indicies_in_agreement);
    points_in_domain = input_points(indicies_of_points,:);


    unit_tangent_vector_of_domain = unit_projection_vectors(index_max_agreement,:);
    unit_orthogonal_vector_of_domain = unit_tangent_vector_of_domain*[0 1; -1 0];
    base_point_of_domain = input_points(combos_paired(index_max_agreement,1),:);

    % Find domain ranges related to point-to-point fit
    projection_vectors_from_base = points_in_domain - base_point_of_domain;
    tangent_distances = sum(projection_vectors_from_base.*unit_tangent_vector_of_domain,2);
    orthogonal_distances = sum(projection_vectors_from_base.*unit_orthogonal_vector_of_domain,2);
   
    % Sort the points in direction of point-to-point projection
    [~,sorted_indicies] = sort(tangent_distances);
    sorted_points_in_domain = points_in_domain(sorted_indicies,:);

    min_tangent_distance = min(tangent_distances); % Can speed this up by using indicies above
    max_tangent_distance = max(tangent_distances); % Can speed this up by using indicies above
    max_orthogonal_distance = max(abs(orthogonal_distances));
    domain_box = ones(5,1)*base_point_of_domain + ...
        [...
        min_tangent_distance*unit_tangent_vector_of_domain - max_orthogonal_distance*unit_orthogonal_vector_of_domain;
        max_tangent_distance*unit_tangent_vector_of_domain - max_orthogonal_distance*unit_orthogonal_vector_of_domain;
        max_tangent_distance*unit_tangent_vector_of_domain + max_orthogonal_distance*unit_orthogonal_vector_of_domain;
        min_tangent_distance*unit_tangent_vector_of_domain + max_orthogonal_distance*unit_orthogonal_vector_of_domain;
        min_tangent_distance*unit_tangent_vector_of_domain - max_orthogonal_distance*unit_orthogonal_vector_of_domain;
        ];

    % Find point-to-point best-fit line
    box_fit_line = ones(2,1)*base_point_of_domain + ...
        [min_tangent_distance*unit_tangent_vector_of_domain; max_tangent_distance*unit_tangent_vector_of_domain];



    % Save domain results
    domains{domain_count}.points_in_domain = points_in_domain; %#ok<SAGROW>
    domains{domain_count}.sorted_points_in_domain = sorted_points_in_domain; %#ok<SAGROW>
    domains{domain_count}.point_to_point_fit_base_point_of_domain = base_point_of_domain; %#ok<SAGROW>
    domains{domain_count}.point_to_point_fit_unit_tangent_vector_of_domain = unit_tangent_vector_of_domain; %#ok<SAGROW>


    % Find regression-fit line
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(sorted_points_in_domain,fig_points_plot);
    
    % Calculalate the domain box for regression fit


    % Plot the results of the point fit
    figure(fig_points_plot);
    h_plot = plot(base_point_of_domain(:,1),base_point_of_domain(:,2),'.','MarkerSize',30);
    current_color = get(h_plot,'Color');
    plot(points_in_domain(:,1),points_in_domain(:,2),'.','MarkerSize',10,'Color',current_color);   
    domainShape = polyshape(domain_box(:,1),domain_box(:,2),'Simplify',false,'KeepCollinearPoints',true);
    plot(domainShape,'FaceColor',current_color,'EdgeColor',current_color,'Linewidth',1);


    % Prune out the points in agreement
    remaining_points(indicies_of_points,:) = nan(size(input_points(indicies_of_points,:)));
    [remaining_agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(remaining_points,unit_projection_vectors, combos_paired, lateral_match_threshold);
    old_remaining_agreements = remaining_agreements;
    remaining_agreements = sum(remaining_agreement_indicies,2);

    % Plot the histogram update
    figure(fig_points_plot+1);
    plot(angles,old_remaining_agreements-remaining_agreements,'.','Color',current_color,'MarkerSize',15);


end

unfitted_points = remaining_points(~isnan(remaining_points(:,1)),:);

% Plot the unfitted agreements
figure(fig_points_plot+1);
nexttile;
plot(angles,remaining_agreements,'k.','MarkerSize',20);
ylim([0 max_original_agreement*1.1]);



% Plot the points that do not fit
figure(fig_points_plot);
plot(unfitted_points(:,1),unfitted_points(:,2),'ro','MarkerSize',10);


%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

%% fcn_INTERNAL_calcUnitVector
function unit_projection = fcn_INTERNAL_calcUnitVector(projection_vector)
vector_length = sum(projection_vector.^2,2).^0.5;
unit_projection = projection_vector./vector_length;
end % Ends fcn_INTERNAL_calcUnitVector


%% fcn_INTERNAL_findAgreementsOfPointsToFits
function [agreement_indicies] = fcn_INTERNAL_findAgreementsOfPointsToFits(input_points,unit_projection_vectors, combos_paired, threshold)

N_combos = length(unit_projection_vectors(:,1));
N_points = length(input_points(:,1));

agreements = zeros(N_combos,1);
agreement_indicies = zeros(N_combos,N_points);

for ith_vector = 1:N_combos
    unit_orthogonal_vector = unit_projection_vectors(ith_vector,:)*[0 1; -1 0];
    base_projection_vectors = input_points - input_points(combos_paired(ith_vector,1),:);
    lateral_distances = sum(unit_orthogonal_vector.*base_projection_vectors,2);

    indicies_in_agreement = (abs(lateral_distances)<threshold)';
    agreement_indicies(ith_vector,:) = indicies_in_agreement;

    agreements(ith_vector) = sum(indicies_in_agreement,2);

end
end % Ends fcn_INTERNAL_findAgreementsOfPointsToFits
