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



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_plotCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***


#### **fcn_geometry_circleCenterFrom3Points**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAngleUsing2PointsOnCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findAngleUsing3PointsOnCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointsFromPointToCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointFromPointToCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_polarLineFrom2PolarCoords**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findVisibleArcsFromPoints**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointsTwoCircles**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findTangentPointTwoCircles**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findPhiConstraints**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findIntersectionLineSegmentWithCircle**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_selfCrossProduct**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_euclideanPointsToPointsDistance**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_findIntersectionOfSegments**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitSlopeInterceptNPoints**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_fitVectorToNPoints**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_flagPointsFurtherFromOriginThanLineSegment**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
for a full test suite.

<a href="#pathplanning_geomtools_geomclasslibrary">Back to top</a>

***

#### **fcn_geometry_flagPointsCloserToOriginThanLineSegment**



**FORMAT:**
```MATLAB
predicted_values = fcn_transform_predictWheelVelocity(pos_rear_left,pos_rear_right,chassis_w,chassis_v)
```

**INPUTS:**

       
**OUTPUTS:**

    
**Dependencies:**


**Examples:**

See the script: script_test_fcn_transform_predictWheelVelocity
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
