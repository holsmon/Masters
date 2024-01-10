%% Alexis Shindhelm 7/31/18; Last updated 8/28/18
%% Hollie Mullin 7/27/2023 Last updated 
% Notes about script:
% 1) Make sure you have loaded in your master_connectome file before
% beginning or that it's in your current folder. Should have one connectome/participant.
% Should be of format connectomeXconnectomeX#participants
% For example: "master connectome" 471X471X136 double means 471X471
% connectome for 136 participants
% 2) Change name of your own specific connectome from 1)
% 3) Make sure you either import correct BCT functions
% or copy to your folder the following: 
% segregation_modified_HM_with_0

 
%% Input csv Files for All Subjects (xcp_d)

%Use code below to feed in matrices. This creates a 3d matrix 
% (ROI x ROI x subject number). This contains each subject's matrix info.

dir_files = dir('/storage/group/nad12/default/fgh3/Users/Hollie/xcp_d/custom_atlases/graph_theory_outputs/power/sess2/rest2/*.csv');
for f=1:length(dir_files) 
    pearson_matrix = readmatrix([dir_files(f).folder filesep dir_files(f).name]);
    % edit out the column headers from .csv file from xcp_d
    pearson_matrix = pearson_matrix(:,:);
    %Perform pearson's r to z transformation. Converts correlation
    %coefficients from xcp_d matrices to z values. 
    r_to_z_matrix = atanh(pearson_matrix);
    PatientTable{1,1,f} = r_to_z_matrix;
end
PatientTable3d = double(cell2mat(PatientTable));

%% Create 3D Matrix of All Subjects with No Negatives (row x column x patient) 
%THIS IS OPTIONAL DEPENDING ON WHETHER YOU WANT NEGATIVES INCLUDED IN YOUR
%DATA!

% We are creating a 3d matrix with negative values set to NaN.  
%0s and negative values are not included in final analysis
MasterMatrix_positive_values = PatientTable3d;

% Make Negative Values Equal to NaN
MasterMatrix_positive_values(MasterMatrix_positive_values < 0)=nan;
%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=MasterMatrix_positive_values;

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas

%% Add in your Atlas 
%This should be 400 regions
%%%MAKE SURE THIS IS SORTED RIGHT based on the order of your correlation
%%%matrices, i.e., master_connectome. 
[Ci,txt,raw] = xlsread('/storage/group/nad12/default/fgh3/Users/Hollie/xcp_d/matlab_graph_theory_output/atlas_labels/Power/Power_atlas_labels_mat.xls');
%inverse of the atlas labels (labels should all be in 1 row, not 1 column -
%left to right)
community_affiliation_vector=Ci';

%% Maximum number of modules/scan
% This should be 17 because based on Yeo atlas and Scafer 400 labels. 
number_of_modules=max(community_affiliation_vector(1,:));

%% Within and Between Connectivity per subject (Chan et al. 2014)
%  Modified to also return subnetworks - MDR (2021)
%This is the correct script for segregation, within, between etc. 

%S, W, and B are the global values for segregation, within network
%connectivity, and between network connectivity

%S_Mod, W_Mod, and B_Mod are provide values for each module for segregation, within, and between
%instead of globally. 

clear s s_nodal S W B S_Mod
for i=1:NP 
    %This is using the segregation_modified script. 
    % This script is calculating the average for within, between,
    % segregation, etc. 
    %Diagonal isn't included in calculation since only upper right portion
    %of matrix triangle is included. 
    %Using MATLAB nanmean in the segregation_modified script since we have
    %NaNs present in our functional connectivity matrices. 
    [S(i,1), W(i,1), B(i,1), S_Mod(i,:), W_Mod(i,:), B_Mod(i,:)] = w_b_s_calc_with_nan(master_connectome(:,:,i),Ci); 
end
