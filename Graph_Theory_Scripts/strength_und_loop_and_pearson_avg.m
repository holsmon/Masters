%% Hollie Mullin 9/16/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

% Note that this script MUST be run WITHOUT negative correlations. See
% strength_und_sign_loop script if you would like to calculate weights
% based on negative and positive values in your matrix. 

 
%% Input csv files for all Subjects

%Use code below to feed in connectivity matrices. This creates a 3d matrix 
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

%% Set diagonal to nan
% for brain connectivity toolbox scripts, diagonals should be set to nan to
% indicate no relationship (1:1 connections between nodes to themselves).

%Otherwise, this script would read the 1 as a meaningful connections
%between the nodes

for n = 1:size(PatientTable3d,1)
    PatientTable3d(n,n,:) = nan;
end

%% Set negative correlations to nan 
% Make sure that negative correlations are NOT connected when evaluating
% strength. 

%Create matrix that represents negatives have been turned into nan
PatientTable3d_positive = PatientTable3d;

% Make Negative Values Equal to nan
PatientTable3d_positive(PatientTable3d_positive < 0) = nan;

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d_positive;% Change master_connectome to connectome matrix

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas

%% Getting strength per subject
%This part of the script calculates strength for every single node. Should
%be an 87 x atlas size (e.g., 400 nodes) output

%formula from strengths_und.m script from Olaf Sporn BCT

%strength_per_node is output variable we care about

for i=1:NP  
    %This is using the strengths_und.m script 
    % This script is calculating average strength of each subject's
    % connectivity matrix 

    %Calculates sum across 1 row/column NOT both (no double calculations
    %here)
    [strength_per_node(i,:)] = nansum(master_connectome(:,:,i)); 
end

%% Average sum of weights across all nodes.

%calculates average of sum of weights for ALL nodes (from
%strength_per_node_variable) for each subject. 

sum_of_nodes=strength_per_node;
%Make 0 (not meaninful/missing data) to nan so they are not affecting the
%average. 
sum_of_nodes(sum_of_nodes <.001) = nan;

%average_strength_across_nodes is output variable we care about here. 

for i=1:NP   
    [average_strength_across_nodes(i,:)] = nanmean(sum_of_nodes(i,:)); 
end


%% Getting average pearson correlation coef per subject
%This script calculates the average pearson correlation for every subject across all ndoes
% . This is averaged across the pearson correlation coeff of all nodes.

% only evaluates the upper triangle of the matrix
%average_correlation is output variable we care about here

for j = 1:NP
    UpperTriangConn=triu(master_connectome(:,:,j));
    average_correlation(j,1)=mean(mean(UpperTriangConn,'omitnan'),'omitnan');
end
