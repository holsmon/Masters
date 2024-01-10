%% Hollie Mullin 9/19/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

% WEIGHT_CONVERSION    Conversion of weights in input matrix
%
%   W_bin = weight_conversion(W, 'binarize');
%   W_nrm = weight_conversion(W, 'normalize');
%   L = weight_conversion(W, 'lengths');
%   W_fix = weight_conversion(W, 'autofix');
%
%   This function may either binarize an input weighted connection matrix,
%   normalize an input weighted connection matrix, convert an input
%   weighted connection matrix to a weighted connection-length matrix, or
%   fix common connection problems in binary or weighted connection matrices.
%
%       Binarization converts all present connection weights to 1.
%
%       Normalization rescales all weight magnitudes to the range [0,1] and
%   should be done prior to computing some weighted measures, such as the
%   weighted clustering coefficient.
%
%       Conversion of connection weights to connection lengths is needed
%   prior to computation of weighted distance-based measures, such as
%   distance and betweenness centrality. In a weighted connection network,
%   higher weights are naturally interpreted as shorter lengths. The
%   connection-lengths matrix here is defined as the inverse of the
%   connection-weights matrix.
%
%       Autofix removes all Inf and NaN values, remove all self connections
%   (sets all weights on the main diagonal to 0), ensures that symmetric matrices
%   are exactly symmetric (by correcting for round-off error), and ensures that
%   binary matrices are exactly binary (by correcting for round-off error).
%
%   Inputs: W           binary or weighted connectivity matrix
%           wcm         weight-conversion command - possible values:
%                           'binarize'      binarize weights
%                           'normalize'     normalize weights

%                           'lengths'       convert weights to lengths
%                           'autofix'       fixes common weights problems
%
%   Output: W_          output connectivity matrix
%
%
%   Mika Rubinov, U Cambridge, 2012

%   Modification History:
%   Sep 2012: Original
%   Jan 2015: Added autofix feature.
%   Jan 2017: Corrected bug in autofix (thanks to Jeff Spielberg)

%% Input csv files for all Subjects

% Use code below to feed in connectivity matrices. This creates a 3d matrix 
% (ROI x ROI x subject number). This contains each subject's matrix 
dir_files = dir('/storage/group/nad12/default/fgh3/Users/Hollie/xcp_d/custom_atlases/graph_theory_outputs/power/sess1/rest1/*.csv');
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

%% Set nans to 0
% need to set nan (missing data) values to 0 for this script to work. Note
% that setting to 0 is OK here (0's won't affect calculations). This is
% because this script only looks at values greater than or less than 0 in
% calculations. 
PatientTable3d_nan_to_0=PatientTable3d;

PatientTable3d_nan_to_0(isnan(PatientTable3d_nan_to_0))=0;

%% Set diagonal to 0
% for brain connectivity toolbox scripts, diagonals should be set to nan to
% indicate no relationship (1:1 connections between nodes to themselves).

%Otherwise, this script would read the 1 as a meaningful connections
%between the nodes

for n = 1:size(PatientTable3d_nan_to_0,1)
    PatientTable3d_nan_to_0(n,n,:) = 0;
end

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d_nan_to_0;% Change master_connectome to connectome matrix

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas

%% Pick if you want to binarize, normalize, get lengths, or autocorrect

for i =1:NP 

%Choose what you want this script to do. See cases below. 
    if ~exist('wcm','var')
        wcm = 3;
    end

    W=master_connectome(:,:,i);

% loop through all subjects 

%Differnet switch cases
    switch wcm
        case 1
            %binarize
            %if value is not equal to 0, then value is converted to 1
            W=double(W~=0);         % binarize
            
            %Can do something like this. This sets all values >0 = 1.
            %Excludes 0s and negative values. 
            %W=double(W>0);

        case 2
            %'normalize'
            W=W./max(abs(W(:)));    % rescale by maximal weight
        case 3
            %'lengths'
            E=find(W);
            W(E)=1./W(E);           % invert weights
        case 4
            %'autofix'

            % clear diagonal
            n = length(W);
            W(1:n+1:end)=0;
        
            % remove Infs and NaNs. Sets those values to 0. 
            idx = isnan(W) | isinf(W);
            if any(any(idx))
                W(idx)=0;
            end
        
            % ensure exact binariness
            U = unique(W);
            if nnz(U) > 1
                idx_0 = abs(W  ) < 1e-10;
                idx_1 = abs(W-1) < 1e-10;
                if all(all(idx_0 | idx_1))
                    W(idx_0)=0;
                    W(idx_1)=1;
                end
            end
        
            % ensure exact symmetry
            if ~isequal(W,W.')
                if max(max(abs(W-W.'))) < 1e-10
                    W=(W+W).'/2;
                end
            end
        otherwise
            error('Unknown weight-conversion command.')
    end
end
