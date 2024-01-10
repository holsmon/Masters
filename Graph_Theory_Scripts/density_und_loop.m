%% Hollie Mullin 9/17/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

%   kden = density_und(CIJ);
%   [kden,N,K] = density_und(CIJ);
%
%   Density is the fraction of present connections to possible connections.
%
%   Input:      CIJ,    undirected (weighted/binary) connection matrix
%
%   Output:     density
%               node_size (changed from N)  =    number of vertices
%               actual_connections (changed from K) =      number of edges
%
%   Notes:  Assumes CIJ is undirected and has no self-connections.
%           Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008

% Modification history:
% 2009-10: K fixed to sum over one half of CIJ [Tony Herdman, SFU]

 
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
%Density (like degree) only cares about number of links to a node. Ignores connection
%weights. 

%Change threshold at the very end (e.g., .1, .2, .3 etc.)
%Anything under the threshold will NOT have a weight of 1 attached to it. 
for i=1:NP  
    master_connectome_binary(:,:,i) = double(master_connectome(:,:,i)>.1);
end

%% Getting density from binarized matrix

%How many vertices (possible connections) exist? Should be equal to atlas
%ROI size. 
for i=1:NP
    node_size(i,:) = size(master_connectome_binary(:,:,i),1);
end

%Getting how many edges actually exist (how many connections based on
%whatever threshold set). Connections are equal to 0 or 1.
%nnz returns nubmer of non-zero elements in the matrix in upper triangle
for i=1:NP
    actual_connections(i,:) = nnz(triu(master_connectome_binary(:,:,i)));
end

%Finally, calculate density. This is the fraction of actual connections
%divided by the total number of possible connections
for i=1:NP
    density(i,:) = actual_connections(i,:)/((node_size(i,:)^2-node_size(i,:))/2);
end

