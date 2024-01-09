function [best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, flag_is_a_circle] = ...
    fcn_geometry_fitHoughArc(points, transverse_tolerance, varargin)
% fcn_geometry_fitHoughArc
%
% This function is just a renamed version of fitHoughCircle 
% 

[best_fitted_parameters, best_fit_source_indicies, best_agreement_indicies, flag_is_a_circle] = ...
    fcn_geometry_fitHoughCircle(points, transverse_tolerance, varargin);