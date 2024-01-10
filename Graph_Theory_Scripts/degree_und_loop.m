%% Hollie Mullin 9/17/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

%   deg = degrees_und(CIJ);
%
%   Node degree is the number of links connected to the node.
%
%   Input:      CIJ,    undirected (binary/weighted) connection matrix
%
%   Output:     deg,    node degree
%
%   Note: Weight information is discarded.
%   Olaf Sporns, Indiana University, 2002/2006/2008

 
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

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d;% Change master_connectome to connectome matrix

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas


%% Binarize matrix 
%Degree only cares about number of links to a node. Ignores connection
%weights. 

%Change threshold at the very end (e.g., .1, .2, .3 etc.)
%Anything under the threshold will NOT have a weight of 1 attached to it. 
for i=1:NP  
    master_connectome_binary(:,:,i) = double(master_connectome(:,:,i)>.1);
end

%% Getting degree from binarized matrix per node:
%Sums up all the degrees (1) in the matrix to get the total number of
%connections for each node. 

%Sum of degrees per node

for i=1:NP  
    [degrees_per_node(i,:)] = nansum(master_connectome_binary(:,:,i)); 
end

%% Getting degree average across nodes
%Sums up all the degrees (1) in the matrix to get the total number of
%connections for each node. 

sum_of_degrees=degrees_per_node;

%Make 0 (not meaninful/missing data) to nan so they are not affecting the
%average. 
sum_of_degrees(sum_of_degrees <.001) = nan;

for i=1:NP  
    [average_degrees(i,:)] = nanmean(sum_of_degrees(i,:)); 
end

