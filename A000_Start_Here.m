%% HEADER

% A000_Start_Here.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
%
% INPUTS: N/A.
% OUTPUTS:N/A.

% Ver report: (do what you have to)
% MATLAB Version: 9.9.0.1592791 (R2020b) Update 5
% MATLAB License Number: 40469100
% Operating System: Microsoft Windows 11 Home Version 10.0 (Build 22000)
% Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% -----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.9         (R2020b)
% Simulink                                              Version 10.2        (R2020b)
% Bioinformatics Toolbox                                Version 4.15        (R2020b)
% Communications Toolbox                                Version 7.4         (R2020b)
% Curve Fitting Toolbox                                 Version 3.5.12      (R2020b)
% DSP System Toolbox                                    Version 9.11        (R2020b)
% Database Toolbox                                      Version 10.0        (R2020b)
% Image Processing Toolbox                              Version 11.2        (R2020b)
% Mapping Toolbox                                       Version 5.0         (R2020b)
% Partial Differential Equation Toolbox                 Version 3.5         (R2020b)
% Signal Processing Toolbox                             Version 8.5         (R2020b)
% Statistics and Machine Learning Toolbox               Version 12.0        (R2020b)
% Symbolic Math Toolbox                                 Version 8.6         (R2020b)

%% Directions for replicating paper results.

% Figure 1:
    % A: This is an illustration.
    % B: Run A001b_Make_Figure_1_B.m
        % This script requires that you have the 'utils' folder from 
        % https://zenodo.org/record/5140528#.YuA103bMJD9
        % on your path, and you'll likely need permission from the Hantman
        % lab as well.
    % C: Run A001c_Make_Figure_1_C.m
        % This data comes from
        % https://www.nature.com/articles/s41598-018-26780-z
        % You should also seek permission from Ioannis Delis.
    % C: Run A001c_Make_Figure_1_D.m
 
% Figure 2:
    