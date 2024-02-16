# PathPlanning_GeomTools_GeomClassLibrary

<!--
The following template is based on:
Best-README-Template
Search for this, and you will find!
>
<!-- PROJECT LOGO -->
<br />
  <h2 align="center"> PathPlanning_GeomTools_GeomClassLibrary
  </h2>

  <pre align="center">
    <img src=".\Images\Geometry.jpeg" alt="Geom Tools Picture" width="600" height="500">
  </pre>

  <p align="center">
  Library of MATLAB functions related to geometric calculations for paths.
    <br />
  </p>
</p>

***

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about">About</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#repo-structure">Repo Structure</a>
      <ul>
        <li><a href="#directories">Top-Level Directories</li>
        <li><a href="#dependencies">Dependencies</li>
      </ul>
    </li>
    <li><a href="#functions">Functions</li>
        <ul>
          <li><a href="#fcn_geometry_checkinputstofunctions">fcn_geometry_checkInputsToFunctions</li>
          <li><a href="#fcn_geometry_plotcircle">fcn_geometry_plotCircle</li>
          <li><a href="#fcn_geometry_plotarc">fcn_geometry_plotArc</li>
          <li><a href="#fcn_geometry_plotsphere">fcn_geometry_plotSphere</li>
          <li><a href="#fcn_geometry_plotfitdomains">fcn_geometry_plotFitDomains</li>
          <li><a href="#fcn_geometry_calcunitvector">fcn_geometry_calcUnitVector</li>
          <li><a href="#fcn_geometry_calcorthogonalvector">fcn_geometry_calcOrthogonalVector</li>
          <li><a href="#fcn_geometry_shufflepointordering">fcn_geometry_shufflePointOrdering</li>
          <li><a href="#fcn_geometry_circlecenterfrom3points">fcn_geometry_circleCenterFrom3Points</li>
          <li><a href="#fcn_geometry_findangleusing2pointsoncircle">fcn_geometry_findAngleUsing2PointsOnCircle</li>
          <li><a href="#fcn_geometry_findangleusing3pointsoncircle">fcn_geometry_findAngleUsing3PointsOnCircle</li>
          <li><a href="#fcn_geometry_findtangentpointsfrompointtocircle">fcn_geometry_findTangentPointsFromPointToCircle</li>
          <li><a href="#fcn_geometry_findtangentpointfrompointtocircle">fcn_geometry_findTangentPointFromPointToCircle</li>
          <li><a href="#fcn_geometry_arcanglefrom3points">fcn_geometry_arcAngleFrom3Points</li>
          <li><a href="#fcn_geometry_arcdirectionfrom3points">fcn_geometry_arcDirectionFrom3Points</li>
          <li><a href="#fcn_geometry_polarlinefrom2polarcoords">fcn_geometry_polarLineFrom2PolarCoords</li>
          <li><a href="#fcn_geometry_findvisiblearcsfrompoints">fcn_geometry_findVisibleArcsFromPoints</li>
          <li><a href="#fcn_geometry_findtangentpointstwocircles">fcn_geometry_findTangentPointsTwoCircles</li>
          <li><a href="#fcn_geometry_findtangentpointtwocircles">fcn_geometry_findTangentPointTwoCircles</li>
          <li><a href="#fcn_geometry_findphiconstraints">fcn_geometry_findPhiConstraints</li>
          <li><a href="#fcn_geometry_findpointsinsequence">fcn_geometry_findPointsInSequence</li>
          <li><a href="#fcn_geometry_findintersectionlinesegmentwithcircle">fcn_geometry_findIntersectionLineSegmentWithCircle</li>
          <li><a href="#fcn_geometry_selfcrossproduct">fcn_geometry_selfCrossProduct</li>
          <li><a href="#fcn_geometry_euclideanpointstopointsdistance">fcn_geometry_euclideanPointsToPointsDistance</li>
          <li><a href="#fcn_geometry_euclideanpointtopointsdistance">fcn_geometry_euclideanPointToPointsDistance</li>
          <li><a href="#fcn_geometry_findintersectionofsegments">fcn_geometry_findIntersectionOfSegments</li>
          <li><a href="#fcn_geometry_fitslopeinterceptnpoints">fcn_geometry_fitSlopeInterceptNPoints</li>
          <li><a href="#fcn_geometry_fitvectortonpoints">fcn_geometry_fitVectorToNPoints</li>
          <li><a href="#fcn_geometry_flagpointsfurtherfromoriginthanlinesegment">fcn_geometry_flagPointsFurtherFromOriginThanLineSegment</li>
          <li><a href="#fcn_geometry_flagpointsclosertooriginthanlinesegment">fcn_geometry_flagPointsCloserToOriginThanLineSegment</li>
          <li><a href="#fcn_geometry_filllinetestpoints">fcn_geometry_fillLineTestPoints</li>
          <li><a href="#fcn_geometry_fillcircletestpoints">fcn_geometry_fillCircleTestPoints</li>
          <li><a href="#fcn_geometry_fillarctestpoints">fcn_geometry_fillArcTestPoints</li>
          <li><a href="#fcn_geometry_fillemptydomainstructure">fcn_geometry_fillEmptyDomainStructure</li>
          <li><a href="#fcn_geometry_fillspheretestpoints">fcn_geometry_fillSphereTestPoints</li>
          <li><a href="#fcn_geometry_corruptpointswithoutliers">fcn_geometry_corruptPointsWithOutliers</li>
          <li><a href="#fcn_geometry_domainboxbytype">fcn_geometry_domainBoxByType</li>
          <li><a href="#fcn_geometry_findagreementsofpointstolinevector">fcn_geometry_findAgreementsOfPointsToLineVector</li>
          <li><a href="#fcn_geometry_findagreementsofpointstocircle">fcn_geometry_findAgreementsOfPointsToCircle</li>
          <li><a href="#fcn_geometry_findagreementsofpointstoarc">fcn_geometry_findAgreementsOfPointsToArc</li>
          <li><a href="#fcn_geometry_findanglebetweenangles">fcn_geometry_findAngleBetweenAngles</li>
          <li><a href="#fcn_geometry_findarcagreementindicies">fcn_geometry_findArcAgreementIndicies</li>
          <li><a href="#fcn_geometry_fitcircleregressionfromhoughfit">fcn_geometry_fitCircleRegressionFromHoughFit</li>
          <li><a href="#fcn_geometry_fitarcregressionfromhoughfit">fcn_geometry_fitArcRegressionFromHoughFit</li>
          <li><a href="#fcn_geometry_fithoughline">fcn_geometry_fitHoughLine</li>
          <li><a href="#fcn_geometry_fithoughcircle">fcn_geometry_fitHoughCircle</li>
          <li><a href="#fcn_geometry_fithougharc">fcn_geometry_fitHoughArc</li>
          <li><a href="#fcn_geometry_fitlinearregressionfromhoughfit">fcn_geometry_fitLinearRegressionFromHoughFit</li>
          <li><a href="#fcn_geometry_houghsegmentation">fcn_geometry_HoughSegmentation</li>
          <li><a href="#fcn_geometry_houghregression">fcn_geometry_HoughRegression</li>
          <li><a href="#fcn_geometry_fitplanelinearregression">fcn_geometry_fitPlaneLinearRegression</li>
          <li><a href="#fcn_geometry_fitspherelsqregression">fcn_geometry_FitSphereLSQRegression</li>
          <li><a href="#fcn_geometry_fitrightcone">fcn_geometry_fitRightCone</li>
        </ul>
      </ul>
    <li><a href="#usage">Usage</a></li>
     <ul>
     <li><a href="#general-usage">General Usage</li>
     </ul>
    <li><a href="#license">License</a></li>
     <li><a href="#major-release-versionse">Major Release Versions</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

***

<!-- ABOUT -->
## About

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

Library of MATLAB functions related to geometric calculations for paths.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo

   ```sh
   git clone https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary.git
   ```
3. Run the main code in the root of the folder (script_demo_GeomTools.m). This will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory

4. Confirm it works! Run script_demo_GeomTools. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->

## Repo Structure

### Directories

The following are the top level directories within the repository:
<ul>
 <li>Functions folder: Contains all functions and their test scripts.</li>
 <li>Utilities: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often other cloned repositories.</li>
 <li>Images folder: Contains images used in the README file.</li>
</ul>

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

The dependencies are automatically installed by running the root master script (script_demo_GeomTools.m).

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- FUNCTION DEFINITIONS -->
## Functions


#### **fcn_geometry_checkInputsToFunctions**

Checks the variable types commonly used in the geometry codes to ensure they are correctly formed.

This function is typically called at the top of most functions. The input is a variable and a string defining the "type" of the variable. This function checks to see that they are compatible. For example, say there 'column_vector' type of variables used in the function that is always a N x 1 array; if someone had a variable called "test_example", they could check that this fit the 'column_vector' type by calling fcn_geometry_checkInputsToFunctions(test_example,'column_vector'). This function would then check that the array was N x 1, and if it was not, it would send out an error warning.

**FORMAT:**

```MATLAB
fcn_geometry_checkInputsToFunctions(variable,variable_type_string,(optional_arguments))
```

**INPUTS:**

variable: the variable to check

variable_type_string: a string representing the variable type to
check. The current strings include:
<ul>

'column_of_numbers' - checks that the input type is N x 1 and
is a number. Optional input: an integer forcing the value
of N, giving an error if the input variable does not have
length N.

'2column_of_numbers' - checks that the input type is N x 2 and
is a number. Optional input: an integer forcing the value
of N, giving an error if the input variable does not have
length N. Another optional input is a rwo vector [A B] where,
if B is greater than A, then the vector must be A or longer.
If B is less than A, then the vector must be A or shorter. If
B = A, then the vector must be length A, and no shorter or
greater.

'2or3column_of_numbers'  - checks that the input type is N x 2
or N x 3 and is a number. Optional input: an integer forcing
the value of N, giving an error if the input variable does not
have length N. Another optional input is a rwo vector [A B]
where, if B is greater than A, then the vector must be A or
longer. If B is less than A, then the vector must be A or
shorter. If B = A, then the vector must be length A, and no
shorter or greater.

</ul>


Note that the variable_type_string is not case sensitive: for
example, 'station' and 'Station' or 'STAtion' all give the same
result.
      
**OUTPUTS:**

No explicit outputs, but produces MATLAB error outputs if conditions
not met, with explanation within the error outputs of the problem.
    
**Dependencies:**

Uses MATLABs dbstack feature to trace dependencies 

**Examples:**

See the script: script_test_fcn_geometry_checkInputsToFunctions
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_plotCircle**

This function plots a circle by creating a vector of angles spaced 0.01 radians apart, and plotting this as a line around the perimeter.

**FORMAT:**
```MATLAB
fcn_geometry_plotCircle(centers,radii,(fig_num))
```

**INPUTS:**

centers: an [N x 2] vector in [x y] of the points of circle centers

radii: a [N x 1] vector of the radii of the circles (to avoid
calculation time)

(OPTIONAL INPUTS)

format:
<ul>
A format string, e.g. 'b-', that dictates the plot style or
A color vector, e.g. [1 0 0.23], that dictates the line color
</ul>
fig_num: a figure number to plot results.
       
**OUTPUTS:**

(none)

**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB
%% BASIC example for multiple circles
fig_num = 2;
figure(fig_num); axis square; grid minor;

centers = [1 3; 2 4];
radii = [2; 3];
fcn_geometry_plotCircle(centers,radii,[],fig_num);
```

<pre align="center">
  <img src=".\Images\plotCircle.jpg" alt="fcn_geometry_plotCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_plotCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_plotArc**

This function plots an arc by creating a vector of angles
spaced a fixed angle apart, and plotting this from the start angle to the
end angle. 

NOTE: If the end angle is numerically larger than the start angle, the
plot will be in the positive direction; otherwise, the plot will be in
the negative direction starting at the start angle, proceeding to the end
angle.

NOTE: The user can specify the fixed angle used to space the plotting
points via an optional input, degree_step.

**FORMAT:**

```MATLAB
    arc_points = fcn_geometry_plotArc(...
    centers,...
    radii,...
    start_angle_in_radians, 
    end_angle_in_radians,
    (degree_step),
    (format),
    (fig_num))
```

**INPUTS:**

centers: an [N x 2] vector in [x y] of the points of circle centers

radii: a [N x 1] vector of the radii of the circles (to avoid
calculation time)

start_angle_in_radians: the starting angle of the arc, in radians

end_angle_in_radians: the starting angle of the arc, in radians

(OPTIONAL INPUTS)

degree_step: how many degrees between plotting points. Default is 1
degree.

format:

<ul>
A format string, e.g. 'b-', that dictates the plot style or
A color vector, e.g. [1 0 0.23], that dictates the line color
</ul>

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
      
**OUTPUTS:**

arc_points: the [x y] coordinates of the arc points. If N
centers and radii are given, with N>1, then arc_points will be a
cell array of points.

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

```MATLAB
%% BASIC example for one arc plotting from 0 to 90, negative
fig_num = 2;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = -90 * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);
```

<pre align="center">
  <img src=".\Images\plotArc.jpg" alt="fcn_geometry_plotArc picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_plotArc
for a full test suite

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_plotSphere**

This function plots a sphere by creating a vector of angles
spaced 0.01 radians apart, and plotting this as a line around the
perimeter.

**FORMAT:**
```MATLAB
    XYZ_points = fcn_geometry_plotSphere(...
    centers,...
    radii,...
    (color_vector);
    (fig_num))
```

**INPUTS:**

centers: an [N x 2] vector in [x y] of the points of sphere centers

radii: a [N x 1] vector of the radii of the spheres (to avoid
calculation time)

(OPTIONAL INPUTS)

color_vector: A color vector, e.g. [1 0 0.23], that dictates the
sphere color. 

fig_num: a figure number to plot results.
       
**OUTPUTS:**

sphere_points: the [x y z] coordinates of the sphere points. If N
centers and radii are given, with N>1, then sphere_points will be a
cell array of points.

**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB
%% BASIC example for one circle
fig_num = 1;
figure(fig_num); clf;

center = [1 3 5];
radius = [2]; %#ok<*NBRAK>
color_vector = [];
fcn_geometry_plotSphere(center, radius, color_vector, fig_num);
```

<pre align="center">
  <img src=".\Images\plotSphere.jpg" alt="fcn_geometry_plotSphere picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_plotSphere
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_plotFitDomains**

Given a set of domains obtained by Hough or regression fitting, plots
these. An optional figure number can be given.

**FORMAT:**

```MATLAB
 fcn_geometry_plotFitDomains(domains, (fig_num))
```

**INPUTS:**

domains: a structure that includes the domains to be plotted. See
fcn_geometry_HoughSegmentation for example outputs

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
      
**OUTPUTS:**

(none)

**Dependencies:**

(none)

**Examples:**

```MATLAB

%% Fill in some test data
% This takes a while - it's generated from the test script for Hough
% Segmentation

clear example_domains

if ~exist('example_domains','var')
    % Advanced example 3: find segments within a chevron
    M = 10; % 40 points per meter

    rng(234)
    sigma = 0.02;

    multi_segment_test_points = [];

    seed_points = [0 0; 10 0];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points]; %#ok<*NASGU>

    seed_points = [0 0; 10 5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [2 0; 3 1.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [4 0; 5 2.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [6 0; 7 3.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [8 0; 9 4.5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];

    seed_points = [10 0; 10 5];
    subtest_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);
    multi_segment_test_points = [multi_segment_test_points; subtest_points];
     
    % Add outliers to corrupt the results
    outliers = [10*randn(100,1) 5*randn(100,1)];
    multi_segment_test_points = [multi_segment_test_points; outliers];


    % Call the segmentation function
    fig_num = 3;
    transverse_tolerance = 0.1; % Units are meters
    station_tolerance = 0.2; % Units are meters
    threshold_max_points = 10;
    input_points = multi_segment_test_points;

    example_domains = fcn_geometry_HoughSegmentation(multi_segment_test_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);
end

%% BASIC test of plotting
fig_num = 1234;
figure(fig_num);
clf;
hold on;
grid on;
axis equal
grid minor;

plot(multi_segment_test_points(:,1),multi_segment_test_points(:,2),'k.','MarkerSize',20);

fcn_geometry_plotFitDomains(example_domains, fig_num);

```

<pre align="center">
  <img src=".\Images\plotFitDomains.jpg" alt="fcn_geometry_plotFitDomains picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_plotFitDomains
for a full test suite

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_calcUnitVector**

This function finds the unit vectors associated with a list of input vectors

**FORMAT:**

```MATLAB
 unit_vectors = fcn_geometry_calcUnitVector(input_vectors, (fig_num))
```

**INPUTS:**

input_vectors: a list of Nxm vector where N is the number of vectors
that should be converted into unit-length vector, and m is the
dimension of the vector (typically 2 or 3).

(OPTIONAL INPUTS)

fig_num: a figure number to plot the results.
      
**OUTPUTS:**

unit_vectors: the unit-length vectors corresponding to each input vector.

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

```MATLAB
%% Test 2: many vectors
fig_num = 2;
input_vectors = randn(10,2); 
unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));
```

<pre align="center">
  <img src=".\Images\calcUnitVector.jpg" alt="fcn_geometry_calcUnitVector picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_calcUnitVector
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_calcOrthogonalVector**

Finds the unit vectors orthogonal to an input list of input vectors. For
2D inputs, the orthogonal vectors are always positive direction. For N-D
vectors, the orthogonal vectors are always orthogonal, but have random
rotations around the axis of the input vectors.

**FORMAT:**

```MATLAB
 unit_orthogonal_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, (fig_num))
```

**INPUTS:**

input_vectors: a list of Nxm vector where N is the number of vectors
that should be converted into unit-length orthogonal vector, and m
is the dimension of the vector (typically 2 or 3).

(OPTIONAL INPUTS)

fig_num: a figure number to plot the results.

      
**OUTPUTS:**

unit_orthogonal_vectors: the unit-length vectors orthogonal to each input vector.

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector

**Examples:**

```MATLAB
%% Test 4: a basic test in 3D, many points
fig_num = 4;
figure(fig_num);
clf;

step = 0.05;

input_vectors = (step:step:1)'.*[3 2 4]; 
seed_points = [];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));

```

<pre align="center">
  <img src=".\Images\calcOrthogonalVector.jpg" alt="fcn_geometry_calcOrthogonalVector picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_calcOrthogonalVector
for a full test suite


<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_shufflePointOrdering**

given N points, with N>=2, this function creates a set of M points per unit distance
between these points randomly distributed with variance sigma.

**FORMAT:**

```MATLAB
corrupted_points = fcn_geometry_shufflePointOrdering(input_points,
(probability_of_corruption), (magnitude_of_corruption), (fig_num));
```

**INPUTS:**

input_points: a Nx2 vector where N is the number of points, but at
least 2.

(Optional Inputs)

fig_num: the figure number to use for plotting
      
**OUTPUTS:**

shuffled_points: a list of test points that are shuffled in ordering

**Dependencies:**

none

**Examples:**

```MATLAB
%% Test 1: a basic test 
fig_num = 1;
figure(fig_num);
clf;

seed_points = [2 3; 4 5; 3 2];
M = 10;
sigma = 0.02;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);

shuffled_points = fcn_geometry_shufflePointOrdering(test_points, (fig_num));
```

<pre align="center">
  <img src=".\Images\shufflePointOrdering.jpg" alt="fcn_geometry_shufflePointOrdering picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

See the script: script_test_fcn_geometry_shufflePointOrdering
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_circleCenterFrom3Points**

This function calculates the center of a circle from three points given as vectors in x and y

**FORMAT:**
```MATLAB
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,(fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is at least 3. If N = 3, a circle will be
fit between these threee points, if N = 4 or more, then one circle
will be fit to the first three points, another cicle to the next
three points, etc.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.

       
**OUTPUTS:**

centers: an [(N-2)x1] vector of the centers of the circles, in [x y]

radii: the radius of each the circles, as an [(N-2)x1] vector
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**


```MATLAB

```

<pre align="center">
  <img src=".\Images\circleCenterFrom3Points.jpg" alt="fcn_geometry_circleCenterFrom3Points picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_circleCenterFrom3Points
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAngleUsing2PointsOnCircle**

This function calculates the angle from the start_points location to the end_points, in the direction of the vector given by is_clockwise.

**FORMAT:**
```MATLAB
[angles] = 
fcn_geometry_findAngleUsing2PointsOnCircle(centers,radii,start_points_on_circle,end_points_on_circle,cross_products,varargin)
```
**INPUTS:**

centers: an [N x 2] vector in [x y] of the points of circle centers

radii: a [N x 1] vector of the radii of the circles (to avoid
calculation time)

start_points_on_circle: an [N x 2] vector in [x y] of the points
where sectors start

end_points_on_circle: an [N x 2] vector in [x y] of the points
where sectors end

cross_products: an [N x 1] vector denoting the cross product
direction to follow from input point to output point

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.

**OUTPUTS:**

angles: an [N x 1] vector of the angles, in radians, between input
points and output points
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findAngleUsing2PointsOnCircle.jpg" alt="fcn_geometry_findAngleUsing2PointsOnCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findAngleUsing2PointsOnCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAngleUsing3PointsOnCircle**

This function calculates the angle from the start_points location to the end_points, in the direction of the vector given by is_clockwise.

**FORMAT:**
```MATLAB
[angles,better_angles,better_angle_range,inpoints_are_closer_to_apex] = fcn_geometry_findAngleUsing3PointsOnCircle(apex_points,centers,start_points_on_circle,end_points_on_circle,radii,incoming_source_points,outgoing_destination_points,varargin)
```

**INPUTS:**

apex_points: an [N x 2] vector of X,Y data for each apex point

centers: an [N x 2] vector in [x y] of the points of circle centers

start_points_on_circle: an [N x 2] vector in [x y] of the points
where sectors start

end_points_on_circle: an [N x 2] vector in [x y] of the points
where sectors end

radii: a [N x 1] vector of the radii of the circles (to avoid
calculation time)

incoming_source_points: an [N x 2] vector in [x y] of the points
where the incoming line is originating from

outgoing_destination_points: an [N x 2] vector in [x y] of the points
where the outgoing line segment is going to.


(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

angles

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_findTangentPointFromPointToCircle
fcn_geometry_findTangentPointsFromPointToCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findAngleUsing3PointsOnCircle.jpg" alt="fcn_geometry_findAngleUsing3PointsOnCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findAngleUsing3PointsOnCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointsFromPointToCircle**

This function finds the two tangent points on a circle given a point.

This function calculates the tangent points on circles defined by centers and radii, where tangents pass through location given by points. This function allows vectorization where centers and points can be a vector [x y] where x and y are columns. The radii must be a column vector of same length as centers. Points are also in [x y] format, of same length as centers. The output

**FORMAT:**
```MATLAB
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(centers,radii,points,(fig_num))
```

**INPUTS:**

centers: an [N x 2] vector of X,Y data for each circle center

radii: a [N x 1] vector of radii for each circle center

points: an [N x 2] vector of X,Y data for each point

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

points_tangent: an [2N x 2] vector of X,Y data for each circle,
arranged such that the top points for the N circles appear first,
then the bottom points for the N circles. NaN is returned for any
points that are within a circle.
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**


```MATLAB

```

<pre align="center">
  <img src=".\Images\findTangentPointsFromPointToCircle.jpg" alt="fcn_geometry_findTangentPointsFromPointToCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findTangentPointsFromPointToCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointFromPointToCircle**

This fucntion finds the ONE tangent point on a circle given a point and a cross product.

This calculates the tangent points on circles defined by centers and
radii, where tangents pass through location given by points, and keeps
only the point that has the same sign as the given cross product. This
function allows vectorization where centers and points can be a vector [x
y] where x and y are columns. The radii must be a column vector of same
length as centers. Points are also in [x y] format, of same length as
centers. 

**FORMAT:**
```MATLAB
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(centers,radii,points,cross_product_goal,(fig_num))
```

**INPUTS:**

centers: an [N x 2] vector of X,Y data for each circle center

radii: a [N x 1] vector of radii for each circle center

points: an [N x 2] vector of X,Y data for each point

cross_product_goal: a [N x 1] vector of cross products for each circle center

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

points_tangent: an [N x 2] vector of X,Y data for each circle.

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_findTangentPointsFromPointToCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findTangentPointFromPointToCircle.jpg" alt="fcn_geometry_findTangentPointFromPointToCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findTangentPointFromPointToCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_arcAngleFrom3Points**

This function calculates the angle between three
points such that the angle starts at point1, passes through point2, and
ends at point3. The angles are presented as positive

**FORMAT:**
```MATLAB
[arc_angle_in_radians_1_to_2, arc_angle_in_radians_1_to_3, circle_centers, radii, start_angles_in_radians]  = fcn_geometry_arcAngleFrom3Points(points1, points2, points3,(fig_num))
```

**INPUTS:**

point1, point2, point3: a Nx2 vectors of point pairings. Note: this
function is vectorized so that multiple points can be entered
simultaneously

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

arc_angle_in_radians_1_to_2: an [Nx1] vector containing the angle
subtending from point1 to point2, passing toward point3

arc_angle_in_radians_1_to_3: an [Nx1] vector containing the angle
subtending from point1 to point3, passing through point2

circle_centers: the centers of the circles of the arcs, as [Nx2]
vector

radii: the radii of the arcs as an [Nx1] vector

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points
fcn_geometry_arcDirectionFrom3Points
fcn_geometry_findAngleUsing2PointsOnCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\arcAngleFrom3Points.jpg" alt="fcn_geometry_arcAngleFrom3Points picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_arcAngleFrom3Points
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_arcDirectionFrom3Points**

This function calculates whether or not 3 points
make an arc in the clockwise or counter-clockwise diretion. The direction
is calculated by performing a cross product between the vectors from
points1 to points3, versus the vectors from points1 to points2. The
funtion returns the sign of the result.

**FORMAT:**
```MATLAB
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(point1, point2, point3,(fig_num))
```

**INPUTS:**

point1, point2, point3: a Nx2 vectors of point pairings. Note: this
function is vectorized so that multiple points can be entered
simultaneously

(OPTIONAL INPUTS)

fig_num: a figure number to plot results

       
**OUTPUTS:**

is_counterClockwise: an [Nx1] vector containing 1 if the
corresponding row of points creates a counter-clockwise arc
(positive) or -1 if it would create a clockwise arc (e.g. negative).
It returns 0 if the direction is undefined (if points are repeated,
colinear, etc.).

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
cross
sign

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\arcDirectionFrom3Points.jpg" alt="fcn_geometry_arcDirectionFrom3Points picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_arcDirectionFrom3Points
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_polarLineFrom2PolarCoords**

This function converts two polar points into the polar form of a line, giving phi and rho values for the line.
See: http://www.nabla.hr/Z_MemoHU-015.htm

**FORMAT:**
```MATLAB
[phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points,(fig_num))
```

**INPUTS:**

points: an [2 x 2] vector of [theta1 r1; theta2 r2] which are the
ends of the line segment to be fit with a line

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

[phi,rho] = the phi and rho term in the polar form of a line:
cos(theta - phi) = rho/r, where rho is the distance of the line from
the origin, and phi is the angle that the perpendicular from the
origin to the line makes with the polar axis (e.g. the x-axis)
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\polarLineFrom2PolarCoords.jpg" alt="fcn_geometry_polarLineFrom2PolarCoords picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_polarLineFrom2PolarCoords
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findVisibleArcsFromPoints**

This function finds the amount of arc visible on a circle from a point external to the circle, returning the arc in radians

**FORMAT:**
```MATLAB
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(centers,radii,points,(fig_num))
```

**INPUTS:**

centers: an [N x 2] vector of X,Y data for each circle center

radii: a [N x 1] vector of radii for each circle center

points: an [N x 2] vector of X,Y data for each point

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.

       
**OUTPUTS:**

visible_arc_angles: a [N x 1] vector of the arc angle, in radians,
for each circle, that is visible from the outside points. NaN is
returned for any points that are within a circle.

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findVisibleArcsFromPoints.jpg" alt="fcn_geometry_findVisibleArcsFromPoints picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findVisibleArcsFromPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointsTwoCircles**

This function finds tangent points from one set of circles to another.

This calculates all the tangent points on circles defined by centers and
radii, where tangents can be either inner tangents or outer tangents. An
optional flag can force one or the other to be used.

This function allows vectorization where centers and radii can be a
vector [x y] and [r] respectively, where x, y, and r are equal-length
columns. The number of starting centers must match the number of ending
centers, as the circles are matched to each other when calculating
outputs.


**FORMAT:**
```MATLAB
[points_tangent_start,points_tangent_end] = fcn_geometry_findTangentPointsTwoCircles(centers_start,centers_end,radii_start,radii_end,(flag_inside_or_out),(voting_points_start,voting_points_end),(fig_num))
```

**INPUTS:**

centers_start: an [N x 2] vector of X,Y data for each circle center
to start the tangent line.

centers_end: an [N x 2] vector of X,Y data for each circle center
to end the tangent line.

radii_start: a [N x 1] vector of radii for each starting circle

radii_end: a [N x 1] vector of radii for each starting circle

(OPTIONAL INPUTS)

flag_inside_or_outside: a scalar flag indicating whether to
calculate all the tangent points, or some subset. The flag can be
one of the following:
<ul>
flag = 0 (default) to calculate both inside and outside
tangents

flag = 1 to calculate outside only, 

flag = -1 to calculate inside tangents only
</ul>
voting_points_start,voting_points_end: each is a [N x 2] vector of
X,Y points which allows user to enter voting points, to keep only
tangents whose start and end are closest to the voting points.

fig_num: a figure number to plot results.

       
**OUTPUTS:**

points_tangent_start: an [N x 2] vector of X,Y data for each circle.

points_tangent_end: an [N x 2] vector of X,Y data for each circle.

inconsistent_votes: shows consistency check for the votes as an
array the same length of votes. Inconsistent votes will be flagged
as 1.

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_findTangentPointsFromPointToCircle
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findTangentPointsTwoCircles.jpg" alt="fcn_geometry_findTangentPointsTwoCircles picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findTangentPointsTwoCircles
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointTwoCircles**

This function finds tangent points from one set of circles to another, returning only
the one set of points tht matches the given cross products

This function allows vectorization where centers and radii can be a
vector [x y] and [r] respectively, where x, y, and r are equal-length
columns. The number of starting centers must match the number of ending
centers, as the circles are matched to each other when calculating
outputs. The length of the cross products column vector must also match.

**FORMAT:**
```MATLAB
[points_tangent_start,points_tangent_end] = fcn_geometry_findTangentPointTwoCircles(centers_start,centers_end,radii_start,radii_end,cross_products_start,cross_products_end,(fig_num))
```

**INPUTS:**

centers_start: an [N x 2] vector of X,Y data for each circle center
to start the tangent line.

centers_end: an [N x 2] vector of X,Y data for each circle center
to end the tangent line.

radii_start: a [N x 1] vector of radii for each starting circle

radii_end: a [N x 1] vector of radii for each starting circle

cross_products_start: a [N x 1] vector of cross products from the
tangent lines (in start-to-end directions) to the circle centers for
the start circles. This is used to specify whether there is an inner
or outer tangent  and which points to keep. The result is only ONE
set of start and end tangent points for each circle: the one that
has the matching cross product.

cross_products_end: a [N x 1] vector of cross products from the
tangent lines (in start-to-end directions) to the circle centers for
the end circles. This is used to specify whether there is an inner
or outer tangent  and which points to keep. The result is only ONE
set of start and end tangent points for each circle: the one that
has the matching cross product.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.

       
**OUTPUTS:**

points_tangent_start: an [N x 2] vector of X,Y data for each circle.

points_tangent_end: an [N x 2] vector of X,Y data for each circle.
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_findTangentPointsFromPointToCircle
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findTangentPointTwoCircles.jpg" alt="fcn_geometry_findTangentPointTwoCircles picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findTangentPointTwoCircles
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findPhiConstraints**

This routine finds the minimum and maximum angle between the apex and
centroid points given the path geometry

**FORMAT:**
```MATLAB
function [phi_start,change_in_phi] = fcn_geometry_findPhiConstraints(p_apex,vertex_1,vertex_2,(fig_num))
```

**INPUTS:**

p_apex: [N x 2] list of apex points in an x,y format

vertex_1: [N x 2] list of first adjacent vertex to be encountered by path

vertex_2: [N x 2] list of second adjacent vertex to be encountered by path

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

phi_min: given in radians, the lower limit for the range of angles
that centroid points can be wrt the apex points 

phi_max: given in radians, the upper limit for the range of angles
that centroid points can be wrt the apex points

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findPhiConstraints.jpg" alt="fcn_geometry_findPhiConstraints picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findPhiConstraints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findPointsInSequence**

This function finds the indicies before and after a base point that are within a tolerance distance sequence of the base point

**FORMAT:**
```MATLAB
sequence_indicies = fcn_geometry_findPointsInSequence(input_distances, base_point_index, station_tolerance, (fig_num))
```

**INPUTS:**

input_distances: a Nx1 vector of distances where N is the number of
distances.

base_point_index: the index of the point that anchors the sequence

station_tolerance: the maximum distance allowed between any points
in the sequence either before or after the base point.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

       
**OUTPUTS:**

sequence_indicies: the indicies that meet the sequence agreement
criteria around the base point, including the index of the base
point

**Dependencies:**

(none)

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findPointsInSequence.jpg" alt="fcn_geometry_findPointsInSequence picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findPointsInSequence
for a full test suite

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findIntersectionLineSegmentWithCircle**

This function evaluates an edge
against a circle to determine whether there will be any intersections
between the two.

**FORMAT:**
```MATLAB
[intAngle,intPoint] = fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R)
```

**INPUTS:**

pa: a 1 x 2 vector containing the (x,y) coordinates of one end of
the line segment
pb: a 1 x 2 vector containing the (x,y) coordinates of one end of
the line segment
pc: a 1 x 2 vector containing the (x,y) coordinates of the center of
the circle
R: a scalar parameter containing the radius of the circle
       
**OUTPUTS:**

intAngle: an N x 1 vector of intersection angles, relative to the
positive x-axis, with length 0, 1, or 2, depending on the
geometry
intPoint: an N x 2 vector of intersection points, with length 0, 1,
or 2, depending on the geometry
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findIntersectionLineSegmentWithCircle.jpg" alt="fcn_geometry_findIntersectionLineSegmentWithCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findIntersectionLineSegmentWithCircle.m
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_selfCrossProduct**

This function finds self cross product of a path, e.g. whether at each internal point
of a path, the path bends next to the left (positive cross product) or to
the right (negative cross product).

**FORMAT:**
```MATLAB
[cross_products] = fcn_geometry_selfCrossProduct(path,(fig_num))
```

**INPUTS:**

path: an [N x 2] vector of X,Y data for the path, with N>=3.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

cross_products: a [(N-2) x 1] vector of cross products calculated by
crossing segment 1 to segment 2, segment 2 to segment 3. If there N
points in the original path, there will be N-1 segments and
therefore N-2 cross product results. The cross products are
organized such that the first cross-product is between segments 1
and 2, the second is between 2 and 3, etc.

err: returns 0 (good) if all the cross-products are well defined, 1
if there was an error, e.g. any corners are ill-defined such as
stright lines, zero lengths, or segments bent back on themselves.
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_findTangentPointsFromPointToCircle
      
**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\selfCrossProduct.jpg" alt="fcn_geometry_selfCrossProduct picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_selfCrossProduct
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_euclideanPointsToPointsDistance**

This function calculates the distance(s) between a vector of points, POINTS1, and another vector of
points, POINTS2.

**FORMAT:**
```MATLAB
[DIST] = fcn_geometry_euclideanPointsToPointsDistance(POINTS1,POINTS2,(fig_num))
```

**INPUTS:**

POINTS1: an Nx2 or Nx3 series of xy or xyz points 
in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]

POINTS2: an Nx2 or Nx3 series of xy or xyz points 
in the form: [x1 y1 z1; x2 y2 z2; ... ; xn yn zn]

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

DIST: an N x  1 vector of distances [d1; d2; ... ; dn], where N is
the number of point sets
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\euclideanPointsToPointsDistance.jpg" alt="fcn_geometry_euclideanPointsToPointsDistance picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_euclideanPointsToPointsDistance
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_euclideanPointToPointsDistance**

This function finds distance(s) from one point to another point(s)

**FORMAT:**

```MATLAB
[dist] = fcn_geometry_euclideanPointToPointsDistance(pt1,pt2, (fig_num))
```
**INPUTS:**

pt1,pt2: [NxM set of points]

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

dist: distances between the point(s)
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\euclideanPointToPointsDistance.jpg" alt="fcn_geometry_euclideanPointToPointsDistance picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_euclideanPointToPointsDistance
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findIntersectionOfSegments**

This function calculates hits between a sensor
projection and a set of walls, both specified by start and end points,
returning the distance, and location of the hit, and the wall number of
each hit.

**FORMAT:**
```MATLAB
[distance,location,wall_that_was_hit] = fcn_geometry_findIntersectionOfSegments(wall_start,wall_end,sensor_vector_start,sensor_vector_end,(flag_search_type),(fig_num))
```

**INPUTS:**

wall_start: an N x 2 vector containing the X,Y points of the
starting points of each "wall".

wall_end: an N x 2 vector containing the X,Y points of the
ending points of each "wall".

sensor_vector_start: a 1 x 2 vector containing the X,Y points of the
sensor's start location

sensor_vector_end: a 1 x 2 vector containing the X,Y points of the
sensor's end location

(OPTIONAL INPUTS)
flag_search_type: an integer specifying the type of search.

<ul>

0: return distance and location of the first point of overlap,
only if the given sensor_vector overlaps the wall (this is the
default)

1: return distane and location if any projection of the sensor
vector, in any direction, hits the wall (in other words, if
there is any intersection). Note that distance returned will
be negative if the nearest intersection is in the opposite
direction of the given sensor vector.

2: returns distances and locations of all hits (not just the
first one as in the options above) as M x 1 and M x 2 vectors
respectively, where the M rows represent all the detected
intersections.
</ul>

fig_num: a figure number to plot results. Turns debugging on.

       
**OUTPUTS:**

distance: a N x 1 scalar representing the distance to the closest
intersection of the sensor with the path. NaN is returned if not
detected.

location: a N x 2 vector of the X,Y location of intersection point

wall_that_was_hit: the segment number of the wall that was hit (1 is
the first segment, 2 is the second, etc)
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\findIntersectionOfSegments.jpg" alt="fcn_geometry_findIntersectionOfSegments picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_findIntersectionOfSegments.m
for a full test suite. 

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitSlopeInterceptNPoints**

This function finds the slope and intercept of a line connecting two points

**FORMAT:**
```MATLAB
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points)
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2. 
       
**OUTPUTS:**

slope: a scalar (1x1) representing the slope connecting the two
points
intercept: a scalar (1x1) representing the y-axis intercept of the
line fit
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**


```MATLAB

```

<pre align="center">
  <img src=".\Images\fitSlopeInterceptNPoints.jpg" alt="fcn_geometry_fitSlopeInterceptNPoints picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitSlopeInterceptNPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitVectorToNPoints**

Finds the vector fitting through a set of N points where N>=2. The
solution performs a polar-form regression of the points, and thus works
for point clusters in any orientation. For more details, see:

"Feature Extraction and Scene Interpretation for Map-Based Navigation
and Map Building" by Kai Oliver Arras, Roland Y. Siegwart


**FORMAT:**
```MATLAB
[vector_root, unit_vector] = fcn_geometry_fitVectorToNPoints(points)
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, length N>=2. 

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

vector_root: the [2 x 1] matrix of the (x,y) point representing the
root of the vector

unit_vector: the [2 x 1] matrix of the (deltax,deltay) length of the
unit vector attached to the root.

**Dependencies:**

fcn_geometry_checkInputsToFunctions
fcn_geometry_plotCircle

**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\fitVectorToNPoints.jpg" alt="fcn_geometry_fitVectorToNPoints picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitVectorToNPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_flagPointsFurtherFromOriginThanLineSegment**

This function generates a column vector that is 1 if associated point is further from
origin (outside) a line segment, 0 otherwise. Points must be strictly
outside to count (e.g. cannot be ON the line segment).

The method used is to fit a line to the line segment to calculate the
slope and intercept. It then applies the same slope to the test points,
and calculates their intercepts. If the test-point intercepts are closer
to the origin than the line segment, then they are "within" the segment.
For vertical line segments, a special test is done to see if the x-axis
coordinate is closer to the orgin than the line segment.

Note: the points must be a tolerance beyond the line to be "further", and
this tolerance is a factor of the numerical precision of the environment,
eps, typically 10*eps where eps = 2.2204e-16.


**FORMAT:**
```MATLAB
[point_flags] = fcn_geometry_flagPointsFurtherFromOriginThanLineSegment(segment_points,test_points,(fig_num))
```

**INPUTS:**

segment_points: a 2x2 vector where the first row is the [x y]
location of the start of the segment, the second row is the [x y]
location of the end of the segment

test_points: a Nx2 vector where N is the number of points that are
being checked (N must be at least 1);

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

point_flags: a vector (Nx1) representing whether associated points
are closer to the origin than the line segment
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions


**Examples:**

```MATLAB

```

<pre align="center">
  <img src=".\Images\flagPointsFurtherFromOrigin.jpg" alt="fcn_geometry_flagPointsFurtherFromOrigin picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_flagPointsFurtherFromOrigin.m
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_flagPointsCloserToOriginThanLineSegment**

This function generates a column vector that is 1 if associated point is closer to
origin (inside) a line segment, 0 otherwise. Points must be strictly
inside to count (e.g. cannot be ON the line segment).

The method used is to fit a line to the line segment to calculate the
slope and intercept. It then applies the same slope to the test points,
and calculates their intercepts. If the test-point intercepts are closer
to the origin than the line segment, then they are "within" the segment.
For vertical line segments, a special test is done to see if the x-axis
coordinate is closer to the orgin than the line segment.


**FORMAT:**
```MATLAB
[point_flags] = fcn_geometry_flagPointsCloserToOriginThanLineSegment(segment_points,test_points,(fig_num))
```

**INPUTS:**

segment_points: a 2x2 vector where the first row is the [x y]
location of the start of the segment, the second row is the [x y]
location of the end of the segment

test_points: a Nx2 vector where N is the number of points that are
being checked (N must be at least 1);

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

point_flags: a vector (Nx1) representing whether associated points
are closer to the origin than the line segment
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

See the script: script_test_fcn_geometry_flagPointsCloserThanLineSegment
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fillLineTestPoints**

given N points, with N>=2, this function creates a set of M points per unit distance between these points randomly distributed with variance sigma.

**FORMAT:**
```MATLAB
[test_points] = fcn_geometry_fillLineTestPoints(seed_points, M, sigma)
```

**INPUTS:**

seed_points: a Nx2 or Nx3 vector where N is the number of points,
but at least 2.

M: the number of test points to generate per unit
distance.

sigma: athe standard deviation in points
       
**OUTPUTS:**

test_points: a list of test points used to test regression fitting

true_base_points, true_projection_vectors, true_distances: the true
values of the fitting
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector
fcn_geometry_calcOrthogonalVector

**Examples:**

See the script: script_test_fcn_geometry_fillLineTestPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fillCircleTestPoints**

given N points, with N>=2, this function creates a set of M points per unit distance
around the perimeter of a circle with points randomly distributed
radially, with variance sigma.

**FORMAT:**

```MATLAB
[test_points] = fcn_geometry_fillCircleTestPoints(seed_points, M, sigma)
```

**INPUTS:**

circle_center: a 1x2 vector denoting the [X Y] location of the
center of the circle

circle_radius: a 1x1 vector denoting the radius of the circle

M: the number of test points to generate per unit
distance around the circumference

sigma: the standard deviation in points in the radial direction

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.  
       
**OUTPUTS:**

test_points: a list of test points used to test regression fitting
    
**Dependencies:**

fcn_geometry_checkInputsToFunctions

**Examples:**

See the script: script_test_fcn_geometry_fillCircleTestPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fillArcTestPoints**

given N seed_points, with N>=3, this function creates test_points that follow arcs
consisting of M points per unit distance, by fitting circles between
consecutive pairs of 3 points. The location of the points on the arcs
follows a radius that is corrupted by perturbations randomly, with
perturbations normally distributed about the nominal radius with variance
sigma.

**FORMAT:**
```MATLAB
[test_points] = fcn_geometry_fillArcTestPoints(seed_points, M, sigma)
```

**INPUTS:**

seed_points: a Nx2 vector where N is the number of points, but at
least 2.

M: the number of test points to generate per unit
distance.

sigma: athe standard deviation in points

       
**OUTPUTS:**

test_points: a list of test points used to test regression fitting

true_circle_centers, true_circle_radii: the true values of the
circle equations used to fill the points

    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points
fcn_geometry_arcDirectionFrom3Points
fcn_geometry_calcUnitVector

**Examples:**

See the script: script_test_fcn_geometry_fillArcTestPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fillEmptyDomainStructure**

This function fills in an empty domain-type structure with the minimum fields. These
fields include:

**FORMAT:**

```MATLAB
emptyDomain = fcn_geometry_fillEmptyDomainStructure((fig_num))
```

**INPUTS:**

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

**OUTPUTS:**

emptyDomain: a domain-type structure containing the following
fields:

emptyDomain.best_fit_type = 'empty';

<ul>

This is a string used to denote the fit type: 

'empty' - no fit yet done

'unfitted' - points left over after a fit is done 

'Hough line' - A Hough fit of line type

'Hough segment' - A Hough fit of line segment type

'Hough circle' - A Hough fit of a circle type

'Hough arc' - A Hough fit of arc type

'Vector regression line fit' - a vector fit of line data,
using regression

'Vector regression segment fit' - a vector fit of a line
segment, using regression

</ul>

emptyDomain.points_in_domain = [nan nan];
<ul>
These are the points in [N x 2] or [N x 3] format that were
used to create the domain, particularly the statistics of the
domain
</ul>

emptyDomain.best_fit_parameters = nan;

<ul>
These are the parameters of the fit. They are defined
in each fit type in the following manner:

'empty' - no fit yet done, so [nan]

'unfitted' - points left over after a fit is done, so [nan]

'Hough line' - 

[unit_projection_vector_x,
unit_projection_vector_y,
base_point_x, 
base_point_y, 
]

'Hough segment' - 

[unit_projection_vector_x,
unit_projection_vector_y,
base_point_x, 
base_point_y, 
station_distance_min,
station_distance_max,
]

'Hough circle' - 
[circleCenter_x.
circleCenter_y,
radius]

'Hough arc' - 

[circleCenter_x.
circleCenter_y,
radius,
start_angle_in_radians, 
end_angle_in_radians,
flag_this_is_a_circle
] 

'Vector regression line fit' - 

[unit_projection_vector_x,
unit_projection_vector_y,
base_point_x, 
base_point_y, 
]

'Vector regression segment fit' -

[unit_projection_vector_x,
unit_projection_vector_y,
base_point_x, 
base_point_y, 
station_distance_min,
station_distance_max,
]

'Regression circle' - 
[circleCenter_x.
circleCenter_y,
radius]

'Regression arc' - 

[circleCenter_x.
circleCenter_y,
radius,
start_angle_in_radians, 
end_angle_in_radians,
flag_this_is_a_circle
] 
</ul>

emptyDomain.best_fit_domain_box = [];
<ul>
These are the points in [N x 2] format that define the outer
boundary of the best-fit domain, generally. For regression
fits, it is the 2-sigma box. For Hough fits, it is the
bounding box used for Hough voting.
</ul>
emptyDomain.best_fit_1_sigma_box = [];
<ul>
These are the points in [N x 2] format that define the outer
boundary of the domain which is the fit, +/- 1 sigma.
</ul>
emptyDomain.best_fit_2_sigma_box = [];
<ul>
These are the points in [N x 2] format that define the outer
boundary of the domain which is the fit, +/- 2 sigma.
</ul>
emptyDomain.best_fit_3_sigma_box = [];
<ul>
These are the points in [N x 2] format that define the outer
boundary of the domain which is the fit, +/- 3 sigma.
</ul>

    
**Dependencies:**

(none)

**Examples:**

See the script: script_test_fcn_geometry_fillEmptyDomainStructure
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fillSphereTestPoints**

given N points, with N>=3, this fuunction creates a set of quasi-uniformly sampled
points around the perimeter of a sphere with points randomly distributed
radially, with variance sigma.


**FORMAT:**
```MATLAB
test_points = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma, (fig_num))
```

**INPUTS:**

N_points: the number of points to generate, with N>=3

sphere_center: a 1x2 vector denoting the [X Y] location of the
center of the sphere

sphere_radius: a 1x1 vector denoting the radius of the sphere

sigma: the standard deviation in points in the radial direction

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.  
       
**OUTPUTS:**

test_points: a list of test points used to test regression fitting
of spheres
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

See the script: script_test_fcn_geometry_fillSphereTestPoints
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_corruptPointsWithOutliers**

given a set of [Nx2] points, this function randomly adds outliers in the orthogonal direction using a random-normal magnitude.

**FORMAT:**
```MATLAB
corrupted_points = fcn_geometry_corruptPointsWithOutliers(input_points,
(probability_of_corruption), (magnitude_of_corruption), (fig_num));
```

**INPUTS:**

input_points: a Nx2 vector where N is the number of points, but at
least 2.

(Optional Inputs)

probability_of_corruption: the probabiity that a given point is an
outlier, from 0 to 1 (default is 0.02)

magnitude_of_corruption: the magnitude of corruption wherein the
outlier multiplied by a random-normal distribution. The default is
2.

fig_num: the figure number to use for plotting
       
**OUTPUTS:**

corrupted_points: a list of test points that are corrupted
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector
fcn_geometry_calcOrthogonalVector

**Examples:**

See the script: script_test_fcn_geometry_corruptPointsWithOutliers
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_domainBoxByType**

This function generates the domain box given a
particular type of domain. NOTE: the inputs change depending on the type
of domain 


**FORMAT:**
```MATLAB
domain_box = fcn_geometry_domainBoxByType(...
type_of_domain,...
(options),
(fig_num))
```

**INPUTS:**

type_of_domain: a string indicating the type of domain to calculate.

Types include:
<ul>

'arc': produces an arc type. If the inputs cause negative radii,
the arc is limited to zero radius. The call is:

domain_box = fcn_geometry_domainBoxByType(...
'arc',...
circleCenter, circleRadius, angles, distance_from_circle_to_boundary,
(fig_num))

</ul>

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
 
       
**OUTPUTS:**

domainShape: a polyshape created from the [x y] coordinates of the
bounding box representing the domain.
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

See the script: script_test_fcn_geometry_domainBoxByType
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAgreementsOfPointsToLineVector**

Given a set of XY points, a unit projection vector and a base point that the vector is
attached to, find the indicies of the points that are within a
transverse_tolerance distance away from the vector. Then, among these
points, finds the indicies that are within a station_tolerance distance
from each other in a group that includes the base point.


**FORMAT:**
```MATLAB
agreement_indicies = fcn_geometry_findAgreementsOfPointsToLineVector( points, unit_projection_vector, base_point_index, transverse_tolerance, station_tolerance, (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2
rows.

unit_projection_vector: a 1x2 vector that is a unit vector in the
direction to check.

base_point_index: the index of the point to use as a base point,
among the input points.

transverse_tolerance: the orthogonal distance between the points and
the linear vector fit that indicate whether a point "belongs" to the
fit. A point belongs to the fit if the transverse distance is less
than or equal to the transverse_tolerance

station_tolerance: the projection distance between the points in a
vector fit, along the direction of the line, that indicate whether a
point "belongs" to the linear fit of a line segment or not. A line
segment is considered continous if station interval distance is less
than or equal to the station_tolerance. Set station_tolerance to an
empty value, [], to avoid line segment fitting and instead find pure
line fits along the vector.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

agreement_indicies: the indicies of the points that are within
agreement of the best-fit parameters, given the transverse and
station tolerance settings. If a station_tolerance is specified, the
indicies are sorted in station-increasing order.

station_distances: the station distances, relative to the base
point, that are in agreement

    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_findPointsInSequence

**Examples:**

See the script: script_test_fcn_geometry_findAgreementsOfPointsToLineVector
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAgreementsOfPointsToCircle**

Given a set of XY points, a circleCenter, and a circleRadius, finds the indicies of the points that are within a
transverse_tolerance distance away from the circle.

**FORMAT:**
```MATLAB
agreement_indicies = fcn_geometry_findAgreementsOfPointsToCircle(points, circleCenter, circleRadius, transverse_tolerance, (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2
rows.

circleCenter: a 1x2 vector that is a the [x y] coordinates of the
circle center to test.

circleRadius: the radius of the circle to test

transverse_tolerance: the orthogonal distance between the points and
the linear vector fit that indicate whether a point "belongs" to the
fit. A point belongs to the fit if the transverse distance is less
than or equal to the transverse_tolerance

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

agreement_indicies: the indicies of the points that are within
agreement of the best-fit parameters, given the transverse and
station tolerance settings. If a station_tolerance is specified, the
indicies are sorted in station-increasing order.
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_findPointsInSequence

**Examples:**

See the script: script_test_fcn_geometry_findAgreementsOfPointsToCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAgreementsOfPointsToArc**

Given a set of XY points, a base_point_index where an arc is rooted, a
circleCenter, and a circleRadius, finds the indicies of the points that
are within a transverse_tolerance distance away from the circle radius,
while keeping within station_tolerance distance from each point within
the cluster centered at the base_point_index.


**FORMAT:**
```MATLAB
agreement_indicies = fcn_geometry_findAgreementsOfPointsToArc(points, circleCenter, circleRadius, transverse_tolerance, (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2
rows.

base_point_index: the index of the point to use as a base point,
among the input points.

circleCenter: a 1x2 vector that is a the [x y] coordinates of the
circle center to test.

circleRadius: the radius of the circle to test

transverse_tolerance: the orthogonal distance between the points and
the linear vector fit that indicate whether a point "belongs" to the
fit. A point belongs to the fit if the transverse distance is less
than or equal to the transverse_tolerance

(OPTIONAL INPUTS)

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the line, that indicate whether a
point "belongs" to the circle fit (if distance is less than or equal
to the tolerance), or is "outside" the fit (if distance is greater
than the tolerance). If left empty, then a circle fit is performed
using only the transverse_tolerance.

flag_force_circle_fit: specify that only circle fits are allowed.
Can be set to:
<ul>
0: search for both arcs and circles (default)

1: search only for circle fits. If station tolerance is given,
then a circle will only be returned if all points meet both the
transverse and station tolerances. Note: to force a circle fit
irregardless of station tolerance, set station_tolerance to an
empty value, e.g. station_tolerance = [];
</ul>

threshold_to_check_arc: a positive integer that represents the
minimum number of points that must be in circle agreement before an
arc is checked. The circle testing is very fast; however, the arc
testing is slow. If the threshold is known beforehand, then entering
this threshold significantly speeds up the analysis.

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**


agreement_indicies: the indicies of the points that are within
agreement of the best-fit parameters, given the transverse and
station tolerance settings. If a station_tolerance is specified, the
indicies are sorted in station-increasing order.

flag_is_a_circle: 1 if the result forms a circle, 0 otherwise

start_angle_in_radians: the angle that the arc starts at, between 0
and 2pi

end_angle_in_radians: the angle that the arc ends at, between 0
and 2pi
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_findAgreementsOfPointsToCircle
fcn_geometry_findArcAgreementIndicies

**Examples:**

See the script: script_test_fcn_geometry_findAgreementsOfPointsToArc
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAngleBetweenAngles**

This function checks whether angles lie between two
different angles

**FORMAT:**
```MATLAB
[isAngleBetween]  = fcn_geometry_findAngleBetweenAngles(start_angle_in_radians, end_angle_in_radians, direction, angles_to_test_in_radians, (fig_num))
```

**INPUTS:**

start_angle_in_radians: the start angle in radians

end_angle_in_radians: the end angle in radians

direction: the direction connecting the start and end angles to
check. Enter 1 for clockwise, anything else for counter-clockwise.

angles_to_test: a vector of [N x 1] angles in radians to check

(OPTIONAL INPUTS)

fig_num: a figure number to plot results.
       
**OUTPUTS:**

isAngleBetween: an [Nx1] vector containing the 1 if the respective test angle
is between the start and end angles, 0 if not
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points
fcn_geometry_arcDirectionFrom3Points
fcn_geometry_findAngleUsing2PointsOnCircle

**Examples:**

See the script: script_test_fcn_geometry_findAngleBetweenAngles
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findArcAgreementIndicies**

Given a set of points, a circle center and radius, and the index of a
source point within the set of points, finds the indicies of the points
that are within a station_tolerance distance from each other, as measured
in arc length along the circle, around the source point. This is useful
to determine sets of points that fall on or around a circle perimeter,
centered near a test point, that form a "contiguous" arc but perhaps not
completely around the circle. Contiguous is defined as having an arc
distance, when the points are projected onto the circle, less than or
equal to the station tolerance.

**FORMAT:**
```MATLAB
[indicies_in_station_agreement, flag_is_a_circle, start_angle_in_radians, end_angle_in_radians] = ...
fcn_geometry_findArcAgreementIndicies(points, circleCenter, circleRadius, index_source_point, station_tolerance, (fig_num));
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2 rows. 

circleCenter: an 1x2 vector containing the [X Y] location of the
circle used to define the arc

circleRadius: an 1x1 scalar defining the radius of the circle used
to define the arc

index_source_point: the index of the point that will be the "root"
point for searching for a contiguous arc among all the points. This
is needed for the common situations where point data may include
many different arc portions along the same circle, and thus one must
define which arc is desired.

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the line, that indicate whether a
point "belongs" to the arc (if distance is less than or equal
to the tolerance), or is "outside" the fit (if distance is greater
than the tolerance). 

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**


indicies_in_station_agreement: the indicies on the points that are
in a contiguous arc that contains the source_point

flag_is_a_circle: a flag that is 1 if the arc extends completely
around the perimeter thus forming a complete circle. It is 0
otherwise.

start_angle_in_radians: the start angle of the arc

end_angle_in_radians: the end angle of the arc.

    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points

**Examples:**

See the script: script_test_fcn_geometry_findArcAgreementIndicies
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitCircleRegressionFromHoughFit**

Given a set of points that are matched to a circle via a Hough vote,
finds the circular regression fit circle and domain box. 

**FORMAT:**
```MATLAB
[regression_fit_circle_center_and_radius, domain_box] = fcn_geometry_fitCircleRegressionFromHoughFit(source_points,associated_points_in_domain, (fig_num))
```

**INPUTS:**

source_points: a 3x2 matrix of the points used to create the Hough circle fit (used to find direction): 

[point1_x  point1_y; point2_x  point2_y; point3_x  point3_y;]

These points are used to determine the direction of the resulting
circle fit.

associated_points_in_domain: an Nx2 list of points that should be
fit with regression, identified as within the domain according to
Hough voting.

(OPTIONAL INPUTS)

fig_num: a figure number to plot the results.
       
**OUTPUTS:**

regression_fit_circle_center_and_radius: a 3x1 matrix of the format:

[circleCenter_x  circleCenter_y; circleRadius]

domain_box: the box that encloses the 2-standard-deviation interval
around the regression circle fit.

radial_errors: the individual errors in each point, radially, in an
[N x 1] matrix

standard_deviation: the standard deviation in the errors

    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector
fcn_geometry_fitSlopeInterceptNPoints

**Examples:**

```MATLAB
%% Fill test data 
fig_num = 21;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

% circle
circle_center = [3 4];
circle_radius = 2;
M = 50; % points per meter
sigma = 0.02;

circle_test_points = fcn_geometry_fillCircleTestPoints(circle_center, circle_radius, M, sigma); % (fig_num));

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_circle_test_points = fcn_geometry_corruptPointsWithOutliers(circle_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

```

<pre align="center">
  <img src=".\Images\fitCircleRegressionFromHoughFit_1.jpg" alt="fcn_geometry_fitCircleRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

```MATLAB
% Basic call with clean data
fig_num = 1;
figure(fig_num);
clf;
hold on;

[regression_fit_circle, domain_box, radial_errors, standard_deviation] = fcn_geometry_fitCircleRegressionFromHoughFit([circle_test_points(1,:); circle_test_points(2,:); circle_test_points(end,:)],circle_test_points, fig_num);

```

<pre align="center">
  <img src=".\Images\fitCircleRegressionFromHoughFit_2.jpg" alt="fcn_geometry_fitCircleRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitCircleRegressionFromHoughFit
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitArcRegressionFromHoughFit**

Given a domain containing a set of points that are matched to an arc via
a Hough vote, finds the arc regression fit and domain box

**FORMAT:**
```MATLAB
[regression_domain, std_dev_transverse_distance] = fcn_geometry_fitArcRegressionFromHoughFit(Hough_domain, (fig_num))
```

**INPUTS:**

Hough_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

regression_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

std_dev_orthogonal_distance: the standard deviation in the point
fit, as measured in the transverse direction (orthogonal to the line
fit). E.g., this is the total-least-squares standard deviation.
    
**Dependencies:**

fcn_geometry_fitCircleRegressionFromHoughFit

**Examples:**


```MATLAB

%% Filling test data for arcs
arc_seed_points = [2 3; 4 5; 6 3];
[arc_true_circleCenter, arc_true_circleRadius] = fcn_geometry_circleCenterFrom3Points(arc_seed_points(1,:),arc_seed_points(2,:),arc_seed_points(3,:),-1);

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(arc_seed_points, M, sigma); %, fig_num);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Fill test data - 2 arcs
twoarc_test_points = [onearc_test_points(1:30,:); onearc_test_points(50:60,:)];
corrupted_twoarc_test_points = [corrupted_onearc_test_points(1:30,:); corrupted_onearc_test_points(50:60,:)];

```

<pre align="center">
  <img src=".\Images\fitArcRegressionFromHoughFit_1.jpg" alt="fcn_geometry_fitArcRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

```MATLAB
%% Fit the onarc_test_points
fig_num = 234;
figure(fig_num); clf;

transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = 20;
flag_force_circle_fit = 0;
expected_radii_range = [1 10];
flag_find_only_best_agreement = [];
flag_use_permutations = [];

fig_num = 235;
figure(fig_num); clf;
inputPoints = corrupted_twoarc_test_points;
domains_corrupted_twoarc_test_points  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));

```
<pre align="center">
  <img src=".\Images\fitArcRegressionFromHoughFit_2.jpg" alt="fcn_geometry_fitArcRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

```MATLAB
%% Basic call with bad data
fig_num = 2;
figure(fig_num);
clf;
hold on;

regression_domain  =  ...
    fcn_geometry_fitArcRegressionFromHoughFit(domains_corrupted_twoarc_test_points{1}, fig_num); 

```

<pre align="center">
  <img src=".\Images\fitArcRegressionFromHoughFit_3.jpg" alt="fcn_geometry_fitArcRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, 
See the script: script_test_fcn_geometry_fitArcRegressionFromHoughFit
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***


#### **fcn_geometry_fitHoughLine**

Checks all permutations between points to fit a line (N choose 2), then
calculates the line fit in polar form (rho and phi), and determines which
of the points are within a tolerance distance of that line fit. This is
calculated for all possible 2-point permutations, ordered in N-choose-2
format. The best-fit parameters are returned

**FORMAT:**
```MATLAB
[best_fitted_parameters, agreement_indicies] = fcn_geometry_fitHoughLine(points,varargin)
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2 rows.

transverse_tolerance: the orthogonal distance between the points and
the linear curve fit that indicate whether a point "belongs" to the
fit. A point belongs to the fit if the transverse distance is less
than or equal to the transverse_tolerance

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the line, that indicate whether a
point "belongs" to the linear fit of a line segment or not. A line
segment is considered continous if station interval distance is less
than or equal to the station_tolerance. Set station_tolerance to an
empty value, [], to avoid line segment fitting and instead find pure
line fits.

(OPTIONAL INPUTS)

points_required_for_agreement: the number of points required for an
agreement to be valid, with minimum value of 3. If left empty, then
the best agreement will always be returned. If a value is given,
line fitting will continue by clustering data until there are no
fits greater than or equal to points_required_for_agreement. The
results of the line fit will be saved in a cell array.

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

domains: a structure that records details of the domain of fitting.
See fcn_geometry_fillEmptyDomainStructure for details.
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector
fcn_geometry_findAgreementsOfPointsToLineVector    
fcn_geometry_fillEmptyDomainStructure
fcn_geometry_plotFitDomains

**Examples:**

```MATLAB
%% Fill in some test data
rng(383);

% Fill test data - 3 segments
seed_points = [2 3; 4 5; 8 0; 9 3]; 
M = 10;
sigma = 0.05;

test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma);

% Add outliers?
% Corrupt the results
probability_of_corruption = 0.1;
magnitude_of_corruption = 4;

test_points = fcn_geometry_corruptPointsWithOutliers(test_points,...
    (probability_of_corruption), (magnitude_of_corruption));


% Shuffle points?
test_points = fcn_geometry_shufflePointOrdering(test_points);

%% Test 2: a basic test of line segment fitting, noisy points
fig_num = 2;
figure(fig_num);
clf;

transverse_tolerance = 0.2;
station_tolerance = 0.4;
points_required_for_agreement = [];

domains= fcn_geometry_fitHoughLine(test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);

```


<pre align="center">
  <img src=".\Images\fitHoughLine.jpg" alt="fcn_geometry_fitHoughLine picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitHoughLine
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitHoughCircle**

This function takes the input points and tolerance as the input and
outputs the fitted parameters and agreement indices of the top-voted
circle

**FORMAT:**
```MATLAB
domains = fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2 rows. 

transverse_tolerance: the orthogonal distance between the points and
the linear curve fit that indicate whether a point "belongs" to the
fit (if distance is less than or equal to the tolerance), or is
"outside" the fit (if distance is greater than the tolerance).

(OPTIONAL INPUTS)

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the line, that indicate whether a
point "belongs" to the circle fit (if distance is less than or equal
to the tolerance), or is "outside" the fit (if distance is greater
than the tolerance). If left empty, then a circle fit is performed
using only the transverse_tolerance.

points_required_for_agreement: the number of points required for an
agreement to be valid, with minimum value of 3. Default is 10. If
left empty, then the best agreement will always be returned. If a
value is given, line fitting will continue by clustering data until
there are no fits greater than or equal to
points_required_for_agreement. The results of the line fit will be
saved in a cell array.

flag_force_circle_fit: specify that only circle fits are allowed.
Can be set to:
<ul>
0:search for both arcs and circles (default)

1:search only for circle fits. If station tolerance is given,
then a circle will only be returned if all points meet both the
transverse and station tolerances. Note: to force a circle fit
irregardless of station tolerance, set station_tolerance to an
empty value, e.g. station_tolerance = [];
</ul>
expected_radii_range: a vector in form of [r_min r_max] indicating
expected radius. Any radii outside this range will not be assessed. 

flag_find_only_best_agreement: set to 1 if want to only keep best
agreement. Otherwise, searches will be continued until none are left
that have more than points_required_for_agreement. Default is 0, to
find all agreements

flag_use_permutations: specify permutation type. Can be set to:
<ul>
0:search for permutations via ordering 1-2-3, then 2-3-4, etc.
This assumes best-fit circle will be formed by points in direct
sequence, and the points are in order (VERY fast)

1: (default) search for permutations via nchoosek, which assumes best-fit
circle will not be formed by points not in direct sequence, but
the points are in order. (somewhat slow, especially in number of
points is greater than 100).  

N: search for permutations by randomly selecting N of them to
test. This essentially turns the Hough transform vote into
RANSAC.
</ul>
fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

       
**OUTPUTS:**

domains: a structure that records details of the domain of fitting.
See fcn_geometry_fillEmptyDomainStructure for details
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points
fcn_geometry_plotCircle 
fcn_geometry_findArcAgreementIndicies

**Examples:**

```MATLAB
%% Fill test data 

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

```

<pre align="center">
  <img src=".\Images\fitHoughCircle_1.jpg" alt="fcn_geometry_fitHoughCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>


```MATLAB
%% BASIC call with arc data, fitting it with a circle by not specifying station tolerance
fig_num = 111;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = [];
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num));


```

<pre align="center">
  <img src=".\Images\fitHoughCircle_2.jpg" alt="fcn_geometry_fitHoughCircle picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitHoughCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitHoughArc**

This function is just a renamed version of fitHoughCircle 

**FORMAT:**
```MATLAB
domains = fcn_geometry_fitHoughCircle(points, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range), (flag_find_only_best_agreement), (flag_use_permutations), (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2 rows. 

transverse_tolerance: the orthogonal distance between the points and
the linear curve fit that indicate whether a point "belongs" to the
fit (if distance is less than or equal to the tolerance), or is
"outside" the fit (if distance is greater than the tolerance).

(OPTIONAL INPUTS)

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the line, that indicate whether a
point "belongs" to the circle fit (if distance is less than or equal
to the tolerance), or is "outside" the fit (if distance is greater
than the tolerance). If left empty, then a circle fit is performed
using only the transverse_tolerance.

points_required_for_agreement: the number of points required for an
agreement to be valid, with minimum value of 3. Default is 10. If
left empty, then the best agreement will always be returned. If a
value is given, line fitting will continue by clustering data until
there are no fits greater than or equal to
points_required_for_agreement. The results of the line fit will be
saved in a cell array.

flag_force_circle_fit: specify that only circle fits are allowed.
Can be set to:
<ul>
0:search for both arcs and circles (default)

1:search only for circle fits. If station tolerance is given,
then a circle will only be returned if all points meet both the
transverse and station tolerances. Note: to force a circle fit
irregardless of station tolerance, set station_tolerance to an
empty value, e.g. station_tolerance = [];
</ul>
expected_radii_range: a vector in form of [r_min r_max] indicating
expected radius. Any radii outside this range will not be assessed. 

flag_find_only_best_agreement: set to 1 if want to only keep best
agreement. Otherwise, searches will be continued until none are left
that have more than points_required_for_agreement. Default is 0, to
find all agreements

flag_use_permutations: specify permutation type. Can be set to:
<ul>
0:search for permutations via ordering 1-2-3, then 2-3-4, etc.
This assumes best-fit circle will be formed by points in direct
sequence, and the points are in order (VERY fast)

1: (default) search for permutations via nchoosek, which assumes best-fit
circle will not be formed by points not in direct sequence, but
the points are in order. (somewhat slow, especially in number of
points is greater than 100).  

N: search for permutations by randomly selecting N of them to
test. This essentially turns the Hough transform vote into
RANSAC.
</ul>
fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

       
**OUTPUTS:**

domains: a structure that records details of the domain of fitting.
See fcn_geometry_fillEmptyDomainStructure for details
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_circleCenterFrom3Points
fcn_geometry_plotCircle 
fcn_geometry_findArcAgreementIndicies

**Examples:**

```MATLAB
%% Fill test data 

% 1 arc
fig_num = 23;
figure(fig_num);
clf;
hold on;
axis equal
grid on;

seed_points = [2 3; 4 5; 6 3];
[true_circleCenter, true_circleRadius] = fcn_geometry_circleCenterFrom3Points(seed_points(1,:),seed_points(2,:),seed_points(3,:),-1);
trueParameters_onearc_test_points = [true_circleCenter true_circleRadius];

M = 10; % Number of points per meter
sigma = 0.02;

onearc_test_points = fcn_geometry_fillArcTestPoints(seed_points, M, sigma); %, fig_num);

% Corrupt the results
probability_of_corruption = 0.3;
magnitude_of_corruption = 1;

corrupted_onearc_test_points = fcn_geometry_corruptPointsWithOutliers(onearc_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

```

<pre align="center">
  <img src=".\Images\fitHoughCircle_1.jpg" alt="fcn_geometry_fitHoughArc picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

```MATLAB
%% BASIC call with arc data, fitting it with an arc by specifying low station tolerance
fig_num = 222;
figure(fig_num); clf;

inputPoints = corrupted_onearc_test_points;
transverse_tolerance = 0.1;
station_tolerance = 1;
points_required_for_agreement = [];
flag_force_circle_fit = [];
expected_radii_range = [];
flag_find_only_best_agreement = []; flag_use_permutations = [];


domains  = ...
fcn_geometry_fitHoughCircle(inputPoints, transverse_tolerance, ...
        (station_tolerance), (points_required_for_agreement), (flag_force_circle_fit), (expected_radii_range),(flag_find_only_best_agreement),(flag_use_permutations), (fig_num));
```

<pre align="center">
  <img src=".\Images\fitHoughArc.jpg" alt="fcn_geometry_fitHoughArc picture" width="500" height="400">
  <figcaption></figcaption>
</pre>


For more examples, see the script: script_test_fcn_geometry_fitHoughCircle
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***


#### **fcn_geometry_fitLinearRegressionFromHoughFit**

Given a domain containting a
set of points that are matched via a Hough vote, finds the regression fit
vector and domain box. 

NOTE: the vector fit is not the same as a least-squares linear
regression, which minimizes sum-of-squares of the VERTICAL errors between
a line fit and the respective points. The method here is similar to
total-least-squares. Namely, the vector is found that approximately
minimizes the sum-of-squares distance between the vector and the
orthogonal projection to each point. Thus, the fit is minimizing
ORTHOGONAL errors. Unlike linear regression, this method works for data
aligned vertically.


**FORMAT:**
```MATLAB
[regression_fit_line_segment, domain_box] = fcn_geometry_fitLinearRegressionFromHoughFit(Hough_domain, (fig_num))
```

**INPUTS:**

Hough_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

regression_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

std_dev_orthogonal_distance: the standard deviation in the point
fit, as measured in the transverse direction (orthogonal to the line
fit). E.g., this is the total-least-squares standard deviation.
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_calcUnitVector
fcn_geometry_fitVectorToNPoints
fcn_geometry_fitSlopeInterceptNPoints 


**Examples:**

```MATLAB
%% Fill test data 
fig_num = 9999;
figure(fig_num); clf;

% Fill in points
seed_points = [2 3; 4 5; 7 0; 9 5; 9 0];
M = 10;
sigma = 0.02;

line_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma, fig_num);


% Corrupt the results with outliers
probability_of_corruption = 0.2;
magnitude_of_corruption = 4; % 4 times the y-range

corrupted_line_test_points = fcn_geometry_corruptPointsWithOutliers(line_test_points,...
    (probability_of_corruption), (magnitude_of_corruption), (fig_num));

% Shuffle points?
shuffled_corrupted_line_test_points = fcn_geometry_shufflePointOrdering(corrupted_line_test_points);

% Demo Hough line fitting

transverse_tolerance = 0.05;
station_tolerance = 2;
points_required_for_agreement = 20;

domains_line_fitting = fcn_geometry_fitHoughLine(shuffled_corrupted_line_test_points, transverse_tolerance, station_tolerance, points_required_for_agreement, fig_num);

```

<pre align="center">
  <img src=".\Images\fitLinearRegressionFromHoughFit_1.jpg" alt="fcn_geometry_fitLinearRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

<pre align="center">
  <img src=".\Images\fitLinearRegressionFromHoughFit_2.jpg" alt="fcn_geometry_fitLinearRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

```MATLAB

%% Basic call - line fitting
fig_num = 1;
figure(fig_num);
clf;
hold on;

regression_domain = fcn_geometry_fitLinearRegressionFromHoughFit(domains_line_fitting{1}, fig_num);

```

<pre align="center">
  <img src=".\Images\fitLinearRegressionFromHoughFit_3.jpg" alt="fcn_geometry_fitLinearRegressionFromHoughFit picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitLinearRegressionFromHoughFit
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_HoughSegmentation**

Given a set of points, attempts to segment the points into line segments,
circles, arcs, etc. by using a Hough transform methodology. 

The method used here is to search through each of the fitting types one
at a time. Each fit must do several things:

Step 1:  use the minimum model order to find all possible permutations of
points that can create a fit. The minimum model order refers to the
minimum number of points needed for the type of fit Given these points,
find all possible index permutations that can form a complete model. The
minimum number of points (M) for each type of model fit is as follows:

<ul>
line or line segment fit - 2 points minimum

circle or arc fit - 3 points minimum

spiral fit - this requires infinite points, and thus only highly
constrained fits would be possible
</ul>

If there are N points, then there are N choose M permutations possible
assuming M is the minimum number of points for the model fit.  NOTE: the
nchoosek funciton in MATLAB produces these combination sequencies
automatically.

Step 2: Given all possible index combinations of points, find the model
fit and the "votes" for that fit. The model fit is simply the geometric
equation created by the M indicies. For a line, this would be the
equation of the line for the 2 points of a particular index combination.
The votes are determined by creating a tolerance around the fit, a
regrion of agreement. Any of the points that are within this region are
considered associated with each fit. The total count of points in
agreement represent the "votes" for that particular fit.

Step 3: Find the highest votes for a fit. The "best" is given by some
point threshold - for example, 10 points, and any fits with 10 points or
more would be considered "good" fits.

Step 4: For the highest votes, save the relevant information such as the
associated points, the model fit, the shape of the domain of the fit.

The above process is repeated from simplest fit types to the most
complex, until no more fits are possible for a given point threshold.
Usually, this starts with the simplest fit type first (lines). After all
lines are fitted, this proceeds to the next most complicated fit
(circles), then the next most complicated (arcs), etc.

NOTE: The fits are ordered from simplest (and fastest) first, then to
more complex. This is necessary not only for speed, but also because of
degeneracy in the fitting process itself. For example, all line segments
are portions of lines, all lines are circles with infinite radius, all
circles are arcs that extend up to 2*pi, all arcs are spirals with a
radial slope equal to zero, etc. In other words, if one starts attempting
to fit points first using advanced and complex geometric representations
- a spiral, for example - before trying simpler ones first such as a line
segment, then the simpler forms would never be found because all the line
segments would be fit by spirals.


**FORMAT:**
```MATLAB
domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, (fig_num))
```

**INPUTS:**

points: a Nx2 vector where N is the number of points, but at least 2 rows. 

threshold_max_points: the number of points below which a match is no
longer performed.

transverse_tolerance: the orthogonal distance between the points and
a curve fit that indicate whether a point "belongs" to the fit
(if distance is less than or equal to the tolerance), or is
"outside" the fit (if distance is greater than the tolerance).

station_tolerance: the projection distance between the points in a
curve fit, along the direction of the fit, that indicate whether a
point "belongs" to the fit (if distance is less than or equal to the
tolerance), or is "outside" the fit (if distance is greater than the
tolerance).

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

domains: the fitted domains of the segments, ordered from the fit
that has the most points included, to the domain that has the least.
    
**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_fitLinearRegressionFromHoughFit 
fcn_geometry_fitCircleRegressionFromHoughFit
fcn_geometry_fitArcRegressionFromHoughFit
fcn_geometry_plotCircle
fcn_geometry_plotArc
fcn_geometry_fitHoughLine
fcn_geometry_fitHoughCircle

**Examples:**

```MATLAB
%% Fill in test points
% Single segment
seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;


single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);

%% Basic example: find one line segment

fig_num = 10; 
figure(fig_num); clf;

transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = corrupted_single_segment_test_points;

domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

```

<pre align="center">
  <img src=".\Images\HoughSegmentation.jpg" alt="fcn_geometry_HoughSegmentation picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_HoughSegmentation
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_HoughRegression**

Given a set of Hough fits of points in regions, performs regression
fitting on each domain.

**FORMAT:**
```MATLAB
domains = fcn_geometry_HoughRegression(Hough_domains, (fig_num))
```

**INPUTS:**

Hough_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.
       
**OUTPUTS:**

regression_domain: a structure that records details of the domain of
fitting. See fcn_geometry_fillEmptyDomainStructure for details.

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions
fcn_geometry_fitLinearRegressionFromHoughFit 
fcn_geometry_fitCircleRegressionFromHoughFit
fcn_geometry_fitCircleRegressionFromHoughFit
fcn_geometry_plotCircle
fcn_geometry_plotArc
fcn_geometry_fitHoughLine
fcn_geometry_fitHoughCircle

**Examples:**

```MATLAB
%% Basic example: find one line segment

% Single segment
seed_points = [5 0; 15 10];
M = 5; % 10 points per meter
sigma = 0.;
rng(3423)

probability_of_corruption = 0.1;
magnitude_of_corruption = 3;


single_segment_test_points = fcn_geometry_fillLineTestPoints(seed_points, M, sigma,-1);
corrupted_single_segment_test_points = fcn_geometry_corruptPointsWithOutliers(single_segment_test_points,...
    (probability_of_corruption), (magnitude_of_corruption),-1);


fig_num = 10; 
figure(fig_num); clf;

transverse_tolerance = 0.05; % Units are meters
station_tolerance = 1; % Units are meters. Usually station tolerance needs to be larger than transverse tolerance, and it needs to be large enough that it can span gaps in corrupted data
threshold_max_points = 10;
input_points = corrupted_single_segment_test_points;

Hough_domains = fcn_geometry_HoughSegmentation(input_points, threshold_max_points, transverse_tolerance, station_tolerance, fig_num);

% Check the regression fit
regression_domains = fcn_geometry_HoughRegression(Hough_domains, fig_num);
fcn_geometry_plotFitDomains(regression_domains, fig_num+1);


```

<pre align="center">
  <img src=".\Images\HoughRegression.jpg" alt="fcn_geometry_HoughRegression picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_HoughRegression
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitPlaneLinearRegression**

This function fits a plane to the points and finds the coefficients, C1, C2, and C3 in the equation

z = C1*x + C2y + C3

**FORMAT:**
```MATLAB
[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = ...
   fcn_geometry_fitPlaneLinearRegression(points,(fig_num))
```

**INPUTS:**

points: a Nx3 vector where N is the number of points, length N>=3. 

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

**OUTPUTS:**

params: the [1 x 3] matrix of the parameters [C1 C2 C3]

standard_deviation_in_z: the standard deviation in the z-error of
the 

z_fit: the model-fit z values

unit_vector: the unit vector in the direction of <A, B, C> for the
equation: 

      A*(x-x0) + B(y-y0) + C(z-z0) = 0

base_point: the location closest to the origin of the plane

standard_deviation_in_plane_orthogonals: the standard deviation in
the point fitting error in the direction of the unit_vector

plane_distances: the distances of each of the N points to the plane,
measured orthogonally from the plane. Returned as an [N x 1 ] vector.

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

```MATLAB
fig_num = 1;
figure(fig_num);
clf;
rng(1823);

true_parameters = [ 0.1 0.2 3]';
points = randn(100,3);
x = points(:,1);
y = points(:,2);
true_z = [x y ones(size(x))]*true_parameters; % Solve for z vertices data

z = true_z;

[parameters, standard_deviation_in_z, z_fit, unit_vector, base_point, standard_deviation_in_plane_orthogonals] = fcn_geometry_fitPlaneLinearRegression([x y z],fig_num);

```

<pre align="center">
  <img src=".\Images\fitPlaneLinearRegression.jpg" alt="fcn_geometry_fitPlaneLinearRegression picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_fitPlaneLinearRegression
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_FitSphereLSQRegression**

Fitting a sphere to points using least squares based on squared differences 
of squared lengths and square radius

**FORMAT:**
```MATLAB
[C,R,E_total,errors] = fcn_LidarPoseEstimation_FitSphereLSQ(XYZ_array)
```

**INPUTS:**

XYZ_array: an Nx3 array contains X,Y and Z coordinates in Velodyne
Lidar Coordinate, with N>=3.

**OUTPUTS:**

C_sphere : a 1x3 vector contains the coordinate of the center of the sphere
R_sphere : the radius of the sphere
E_total  : the sum of the squared differences of squared lengths 
            and square radius
errors   : the differences of point distances to center to fitted
radius

**Dependencies:**

(none)

**Examples:**

```MATLAB
fig_num = 1;
figure(fig_num); clf;

% Fill in the true values
sphere_center = [6 2 0];
sphere_radius = 1;

% Fill in the test data
N_points = 1000;
sigma = 0.02;
XYZ_array = fcn_geometry_fillSphereTestPoints(N_points, sphere_center, sphere_radius, sigma,-1);
assert(length(XYZ_array(:,1))==N_points);


% Call the fitting function
[C_sphere,R_sphere,~, errors] = fcn_geometry_FitSphereLSQRegression(XYZ_array, fig_num);

```

<pre align="center">
  <img src=".\Images\FitSphereLSQRegression.jpg" alt="fcn_geometry_fFitSphereLSQRegression picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_FitSphereLSQRegression
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***
#### **fcn_geometry_fitRightCone**

This function calculates the ratio of radius to height of the
cone

**FORMAT:**
```MATLAB
cone_parameters = fcn_geometry_fitRightCone(inputPoints,ring_id)
```

**INPUTS:**

inputPoints: a Nx3 vectors of point pairings in the one single scan
layer

ring_id : an integer number indicates the id of the layer, ranging 
from 0 to 15

(OPTIONAL INPUTS)

fig_num: a figure number to plot results. If set to -1, skips any
input checking or debugging, no figures will be generated, and sets
up code to maximize speed.

**OUTPUTS:**

cone_parameters: an [1x2] vector containing ratio of radius to 
height of the cone and the id of the scan layer

**Dependencies:**

fcn_DebugTools_checkInputsToFunctions

**Examples:**

```MATLAB
%% Load Data: Use the rawdata_lidar.mat in the Data folder to start, 
% the loaded data is a struct containing two fields
% rawdata.GPS_SparkFun_Temp_ENU and rawdata.Lidar_pointcloud_cell
data_folder = fullfile("..","Data");
addpath(data_folder)
load('PointCloud_Separated_Data_newTarget.mat')

%% Test 1: a basic test, change the index of scan and ring label if needed
idx_scan = 1;
test_ring = ptCloud_pts_layers_separated_cell{idx_scan}.Ring0;
inputPoints = test_ring(:,1:3);
ring_id = test_ring(1,5);
fig_num = 1;
[cone_parameters,fittedPoints,fitting_result] = fcn_geometry_fitRightCone(inputPoints,ring_id,fig_num);


```

<pre align="center">
  <img src=".\Images\fitRightCone.jpg" alt="fcn_geometry_fitRightCone picture" width="500" height="400">
  <figcaption></figcaption>
</pre>

For more examples, see the script: script_test_fcn_geometry_FitSphereLSQRegression
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

### General Usage

Each of the functions has an associated test script, using the convention

```sh
script_test_fcn_fcnname
```

where fcnname is the function name as listed above.

As well, each of the functions includes a well-documented header that explains inputs and outputs. These are supported by MATLAB's help style so that one can type:

```sh
help fcn_fcnname
```

for any function to view function details.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***
<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

## Major Release Versions

This code is still in development (alpha testing)

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary](https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary)

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
