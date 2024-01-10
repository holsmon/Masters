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
% community_louvain.m, 

%% Input csv files for all Subjects

% Use code below to feed in connectivity matrices. This creates a 3d matrix 
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

%% Set diagonal to 0
% for brain connectivity toolbox scripts, diagonals should be set to nan to
% indicate no relationship (1:1 connections between nodes to themselves).

%Otherwise, this script would read the 1 as a meaningful connections
%between the nodes

for n = 1:size(PatientTable3d,1)
    PatientTable3d(n,n,:) = 0;
end

%% Set negative correlations to 0
% Make sure that negative correlations are NOT connected when evaluating
% strength. 

%Create matrix that represents negatives have been turned into positives
PatientTable3d_positive = PatientTable3d;

% Make Negative Values Equal to nan
PatientTable3d_positive(PatientTable3d_positive < 0) = 0;

%% Set nans to 0
% need to set nan (missing data) values to 0 for this script to work. Note
% that setting to 0 is OK here (0's won't affect calculations). 
PatientTable3d_nan_to_0=PatientTable3d_positive;

PatientTable3d_nan_to_0(isnan(PatientTable3d_nan_to_0))=0;

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d_nan_to_0;% Change master_connectome to connectome matrix

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
% This should be 7 because we're excluding the cerebellum and subcortical. 
number_of_modules=max(community_affiliation_vector(1,:));

%% Make master  community affiliation vector for everyone
for j=1:NP
    master_community_affiliation_matrix(j,:)=community_affiliation_vector;
end

% Calculate Modularity - One value per participant
for j = 1:NP
    %for i = 1:150 % Runs this 150X
        [M, Q] = community_louvain_calcmod_only(master_connectome(:,:,j),1,community_affiliation_vector);
        %inputs are 1) connection matrix, 2) gamma value, and 3) community
        %affiliation vector
        %Q_temp_OA(i,1) = Q; M_temp_OA(:,i) = M;
        %clear M Q
    %end
    %modularity(j,1)=max(Q_temp_OA);
    modularity(j,1)=Q;
    s=strcat('Modularity #', num2str(j)); disp(s)
end
clear i j Q_temp_OA M_temp_OA ind