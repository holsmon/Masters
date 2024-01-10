%% Hollie Mullin 9/16/2023 Last updated 
% Notes about script:
%
% See notes within that script for more details regarding the connectivity
% matrix.

%   C = clustering_coef_wu(W);
%
%   The weighted clustering coefficient is the average "intensity"
%   (geometric mean) of all triangles associated with each node.
%
%   Input:      W,      weighted undirected connection matrix
%                       (all weights must be between 0 and 1)
                        %(NORMALIZE THIS!)
%
%   Output:     C,      clustering coefficient vector
%
%   Note:   All weights must be between 0 and 1.
%           This may be achieved using the weight_conversion.m function,
%           W_nrm = weight_conversion(W, 'normalize');
%
%   Reference: Onnela et al. (2005) Phys Rev E 71:065103
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2015

%   Modification history:
%   2007: original
%   2015: expanded documentation

 
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

%% Set nans to 0
% need to set nan (missing data) values to 0 for this script to work. Note
% that setting to 0 is OK here (0's won't affect calculations). This is
% because this script only looks at values greater than or less than 0 in
% calculations. 
PatientTable3d_nan_to_0=PatientTable3d;

PatientTable3d_nan_to_0(isnan(PatientTable3d_nan_to_0))=0;

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d_nan_to_0;% Change master_connectome to connectome matrix

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas

%% Getting clustering coeff per subject
%This is the Onnela et al., formula, which is the BCT's default. 
%formula from clustering_coef_wu.m script from Olaf Sporn BCT

% We care about:
% C_pos (cluster coef per node, positive values only)
% C_neg (cluster coef per node, negative values only)
% Ctot_pos (mean cluster coef across nodes, positive values only)
% Ctot_neg (mean cluster coef across nodes, negative values only)
%Note that negative values are transformed to absolute value for coef_type
%1 and 2, but not 3. See more details on coef_type below.

for participant=1:NP  
    
    %change the coef_type number if you would like to calculate other
    %correlation coeficients.

%           Desired type of clustering coefficient.
%           Options:  
%           1,  (default) Onnela et al. formula, used in original
%               clustering_coef_wu.m. Computed separately for positive &
%               negative weights.
%           2,  Zhang & Horvath formula, similar to Onnela formula except
%               denominator of Onnela formula relies on binarizing the
%               network whereas this denominator is based on weight value,
%               which reduces the sensitivity of this measure to the
%               weights directly connected to the node of interest.
%               Computed separately for positive & negative weights.
%           3,  Constantini & Perugini's generalization of the Zhang &
%               Horvath formula. This formula takes both positive &
%               negative weights into account simultaneously, & is
%               particularly sensitive to non-redundancy in path
%               information based on sign (i.e., when two weights are
%               positive & one negative, or all three are negative, both of
%               which indicate that the weight of the third path is not
%               redundant information). Produces only one value.
%
    if ~exist('coef_type','var')
    coef_type = 1;
    end

    % Setting W to each connectivity per subject. i = subject number, loops
    % through the 3d master_connectome matrix. 
    W=master_connectome(:,:,participant);
    
    %set diagonal to 0
    n = length(W); 
    W(1:n+1:end) = 0;
    
    %Normalize data
    W=W./max(abs(W(:)));

    n            = length(W);                   %number of nodes

    switch coef_type
        case 1     

            %only looks at positive values in matrix, all others are set to
            %0
            W_pos                = W.*(W>0);

        
            
            % this is a logical argument that sums the total number (e.g., 1) of
            % connections that exist for positive values only.
            % Creates a logicla matrix with only values of 0 or 1 (binary)
            % compared to weighted pearson correlation coef. 

            % Note that sum(W_pos,2) would get you sum of correlation coef
            % for a given node
            K_pos                = sum(W_pos~=0,2);
            
            %what is going on here???? .^ means elemental power...
            cyc3_pos             = diag((W_pos.^(1/3))^3);
            
            %same output as K_pos above, except 0's are set to inf
            K_pos(cyc3_pos == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            
            [C_pos(participant,:)]                = cyc3_pos./(K_pos.*(K_pos-1));         %clustering coefficient
            
            %Set C_pos 0's to nan so that they are not calculated as part
            %of mean
            C_pos_0_to_nan = C_pos;
            C_pos_0_to_nan(C_pos_0_to_nan==0)=nan;
            
            [Ctot_pos(participant,:)]             = nanmean(C_pos_0_to_nan(participant,:));
            
            %Take absolute value of negative correlatioins and only look at
            %negative correlations. 
            W_neg                = -W.*(W<0);
            K_neg                = sum(W_neg~=0,2);
            cyc3_neg             = diag((W_neg.^(1/3))^3);
            K_neg(cyc3_neg == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            [C_neg(participant,:)]                = cyc3_neg./(K_neg.*(K_neg-1));         %clustering coefficient
            
            %Set C_neg 0's to nan so that they are not calculated as part
            %of mean
            C_neg_0_to_nan = C_neg;
            C_neg_0_to_nan(C_neg_0_to_nan==0)=nan;
            
            [Ctot_neg(participant,:)]            = nanmean(C_neg_0_to_nan(participant,:));
        case 2
            W_pos    = W.*(W>0);
            cyc3_pos = zeros(n,1);
            cyc2_pos = zeros(n,1);
            for i = 1:n
                for j = 1:n
                    for q = 1:n
                        cyc3_pos(i) = cyc3_pos(i)+(W_pos(j,i)*W_pos(i,q)*W_pos(j,q));
                        if j~=q
                            cyc2_pos(i) = cyc2_pos(i)+(W_pos(j,i)*W_pos(i,q));
                        end
                    end
                end
            end
            cyc2_pos(cyc3_pos == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            [C_pos(participant,:)]                   = cyc3_pos./cyc2_pos;         %clustering coefficient
            
            %Get C_pos 0 to nan
            C_pos_0_to_nan = C_pos;
            C_pos_0_to_nan(C_pos_0_to_nan==0)=nan;

            [Ctot_pos(participant,:)]                = nanmean(C_pos_0_to_nan(participant,:));
        
            W_neg    = -W.*(W<0);
            cyc3_neg = zeros(n,1);
            cyc2_neg = zeros(n,1);
            for i = 1:n
                for j = 1:n
                    for q = 1:n
                        cyc3_neg(i) = cyc3_neg(i)+(W_neg(j,i)*W_neg(i,q)*W_neg(j,q));
                        if j~=q
                        cyc2_neg(i) = cyc2_neg(i)+(W_neg(j,i)*W_neg(i,q));
                        end
                    end
                end
            end
            cyc2_neg(cyc3_neg == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            C_neg(participant,:)                   = cyc3_neg./cyc2_neg;         %clustering coefficient
           
            %Get C_neg 0 to nan
            C_neg_0_to_nan = C_neg;
            C_neg_0_to_nan(C_neg_0_to_nan==0)=nan;
            
            [Ctot_neg(participant,:)]                = nanmean(C_neg_0_to_nan(participant,:));
        case 3
            cyc3         = zeros(n,1);
            cyc2         = zeros(n,1);
        
            for i = 1:n
                for j = 1:n
                    for q = 1:n
                        cyc3(i) = cyc3(i)+(W(j,i)*W(i,q)*W(j,q));
                        if j~=q
                            cyc2(i) = cyc2(i)+abs(W(j,i)*W(i,q));
                        end
                    end
                end
            end
        
            cyc2(cyc3 == 0) = inf;             %if no 3-cycles exist, make C=0 (via K=inf)
            [C_pos(participant,:)]           = cyc3./cyc2;         %clustering coefficient
            
            %Get C_pos 0 to nan
            C_pos_0_to_nan = C_pos;
            C_pos_0_to_nan(C_pos_0_to_nan==0)=nan;
            
            [Ctot_pos(participant,:)]        = nanmean(C_pos_0_to_nan(participant,:));

            
            [C_neg(participant,:)]           = nanmean(size(C_pos_0_to_nan(participant,:)));

            %Get C_neg 0 to nan
            C_neg_0_to_nan = C_neg;
            C_neg_0_to_nan(C_neg_0_to_nan==0)=nan;

            [Ctot_neg(participant,:)]        = nan(size(Ctot_pos(participant,:)));
    end 
end
