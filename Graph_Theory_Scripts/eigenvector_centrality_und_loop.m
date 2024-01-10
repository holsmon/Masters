%% Hollie Mullin 9/17/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

%   v = eigenvector_centrality_und(CIJ)
%
%   Eigenector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector 
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   Inputs:     CIJ,        binary/weighted undirected adjacency matrix.
%
%   Outputs:      v (eigenvector),        eigenvector associated with the largest
%                           eigenvalue of the adjacency matrix CIJ.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%
%   Contributors:
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012
%   Mika Rubinov, University of Cambridge, 2015

%   MODIFICATION HISTORY
%   2010/2012: original (XNZ, RB)
%   2015: ensure the use of leading eigenvector (MR)

% Low key just watch this YT video because I had no clue what eigenvector
% centrality was until this - Note from Hollie Mullin:
%https://www.youtube.com/watch?v=1S1mD0l9FwU

 
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

%% Calculating eigenvector centrality

%Number of nodes/ROIs
for i=1:NP
    n(i,:) = length(master_connectome(:,:,i));
    if n(i,:) < 1000
        %The columns of V present eigenvectors of A. The main diagonal of D present eigenvalues of A.
        %Eigenvector is direction whereas eigenvalue is strength 
        [V(:,:,i),D(:,:,i)] = eig(master_connectome(:,:,i));
    else
        [V(:,:,i),D(:,:,i)] = eigs(sparse(master_connectome(:,:,i)));
    end
    
    %idx is n (number of nodes)
    [~,idx] = max(diag(D(:,:,i)));


    ec = abs(V(:,idx,i));
    v(i,:) = reshape(ec, length(ec), 1);
end

%% Getting eigenvector average across nodes

sum_of_v=v;

%Make 0 (not meaninful/missing data) to nan so they are not affecting the
%average. 
sum_of_v(sum_of_v <.001) = nan;

for i=1:NP  
    [average_v(i,:)] = nanmean(sum_of_v(i,:)); 
end


