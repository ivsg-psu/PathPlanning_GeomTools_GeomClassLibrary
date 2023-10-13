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
    <img src=".\Images\" alt="Geom Tools Picture" width="577" height="246.67">
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
    <li><a href="structure">Repo Structure</a>
      <ul>
        <li><a href="#directories">Top-Level Directories</li>
        <li><a href="#dependencies">Dependencies</li>
      </ul>
    </li>
    <li><a href="#functions">Functions</li>
        <ul>
          <li><a href="#fcn_geometry_checkinputstofunctions">fcn_geometry_checkInputsToFunctions</li>
          <li><a href="#fcn_geometry_plotcircle">fcn_geometry_plotCircle</li>
          <li><a href="#fcn_geometry_circlecenterfrom3points">fcn_geometry_circleCenterFrom3Points</li>
          <li><a href="#fcn_geometry_findangleusing2pointsoncircle">fcn_geometry_findAngleUsing2PointsOnCircle</li>
          <li><a href="#fcn_geometry_findangleusing3pointsoncircle">fcn_geometry_findAngleUsing3PointsOnCircle</li>
          <li><a href="#fcn_geometry_findtangentpointsfrompointtocircle">fcn_geometry_findTangentPointsFromPointToCircle</li>
          <li><a href="#fcn_geometry_findtangentpointfrompointtocircle">fcn_geometry_findTangentPointFromPointToCircle</li>
          <li><a href="#fcn_geometry_polarlinefrom2polarcoords">fcn_geometry_polarLineFrom2PolarCoords</li>
          <li><a href="#fcn_geometry_findvisiblearcsfrompoints">fcn_geometry_findVisibleArcsFromPoints</li>
          <li><a href="#fcn_geometry_findtangentpointstwocircles">fcn_geometry_findTangentPointsTwoCircles</li>
          <li><a href="#fcn_geometry_findtangentpointtwocircles">fcn_geometry_findTangentPointTwoCircles</li>
          <li><a href="#fcn_geometry_findphiconstraints">fcn_geometry_findPhiConstraints</li>
          <li><a href="#fcn_geometry_findintersectionlinesegmentwithcircle">fcn_geometry_findIntersectionLineSegmentWithCircle</li>
          <li><a href="#fcn_geometry_selfcrossproduct">fcn_geometry_selfCrossProduct</li>
          <li><a href="#fcn_geometry_euclideanpointstopointsdistance">fcn_geometry_euclideanPointsToPointsDistance</li>
          <li><a href="#fcn_geometry_findintersectionofsegments">fcn_geometry_findIntersectionOfSegments</li>
          <li><a href="#fcn_geometry_fitslopeinterceptnpoints">fcn_geometry_fitSlopeInterceptNPoints</li>
          <li><a href="#fcn_geometry_fitvectortonpoints">fcn_geometry_fitVectorToNPoints</li>
          <li><a href="#fcn_geometry_flagpointsfurtherfromoriginthanlinesegment">fcn_geometry_flagPointsFurtherFromOriginThanLineSegment</li>
          <li><a href="#fcn_geometry_flagpointsclosertooriginthanlinesegment">fcn_geometry_flagPointsCloserToOriginThanLineSegment</li>
        </ul>
      </ul>
    <li><a href="#usage">Usage</a></li>
     <ul>
     <li><a href="#general-usage">General Usage</li>
     </ul>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

***

<!-- ABOUT -->
## About

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

Add ABOUT here

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

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

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->
### Directories

The following are the top level directories within the repository:
<ul>
 <li>Functions folder: Contains all functions and their test scripts.</li>
 <li>Utilities: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often other cloned repositories.</li>
 <li>Images folder: Contains images used in the README file.</li>
</ul>

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

The dependencies are automatically installed by running the root master script (script_demo_GeomTools.m).

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

See the script: script_test_fcn_geometry_plotCircle
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

See the script: script_test_fcn_geometry_circleCenterFrom3Points
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

See the script: script_test_fcn_geometry_findAngleUsing2PointsOnCircle
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

See the script: script_test_fcn_geometry_findAngleUsing3PointsOnCircle
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

See the script: script_test_fcn_geometry_findTangentPointsFromPointToCircle
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

See the script: script_test_fcn_geometry_findTangentPointFromPointToCircle
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

See the script: script_test_fcn_geometry_polarLineFrom2PolarCoords
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

See the script: script_test_fcn_geometry_findVisibleArcsFromPoints
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

See the script: script_test_fcn_geometry_findTangentPointsTwoCircles
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

See the script: script_test_fcn_geometry_findTangentPointTwoCircles
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

See the script: script_test_fcn_geometry_findPhiConstraints
for a full test suite.

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

See the script: geometry_findIntersectionLineSegmentWithCircle.m
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

See the script: script_test_fcn_geometry_selfCrossProduct
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

See the script: script_test_fcn_geometry_euclideanPointsToPointsDistance
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

See the script: script_test_fcn_geometry_findIntersectionOfSegments.m
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

See the script: script_test_fcn_geometry_fitSlopeInterceptNPoints
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

See the script: script_test_fcn_geometry_fitVectorToNPoints
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

See the script: script_test_fcn_geometry_flagPointsFurtherFromOrigin.m
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

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***
<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***

## Major Release Versions

This code is still in development (alpha testing)

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***

<!-- CONTACT -->
## Contact

Sean Brennan - sbrennan@psu.edu

Project Link: [https://github.com/ivsg-psu/FeatureExtraction_DataTransforms_TransformClassLibrary](https://github.com/ivsg-psu/FeatureExtraction_DataTransforms_TransformClassLibrary)

<a href="#featureextraction_datatransforms_transformclasslibrary">Back to top</a>

***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
