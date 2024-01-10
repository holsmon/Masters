%% Hollie Mullin 9/20/2023 Last updated 
% Notes about script:

% Must get lengths matrix first from weight_conversion_loop script. Then,
% get distance matrix (calculated from lengths matrix). Then, feed in
% distance matrix to finally obtain path length!

%% MUST USE DISTANCE_WEI BEFORE USING CHARPATH SCRIPT.
%Obtains distance matrices. 

% DISTANCE_WEI       Distance matrix (Dijkstra's algorithm)
%
%   D = distance_wei(L);
%   [D,B] = distance_wei(L);
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path 
%   from node u to node v. The average shortest path length is the 
%   characteristic path length of the network.
%
%   Input:      L,      Directed/undirected connection-length matrix.
%   *** NB: The length matrix L isn't the weights matrix W (see below) ***
%
%   Output:     D,      distance (shortest weighted path) matrix
%               B,      number of edges in shortest weighted path matrix
%
%   Notes:
%       The input matrix must be a connection-length matrix, typically
%   obtained via a mapping from weight to length. For instance, in a
%   weighted correlation network higher correlations are more naturally
%   interpreted as shorter distances and the input matrix should
%   consequently be some inverse of the connectivity matrix. 
%       The number of edges in shortest weighted paths may in general 
%   exceed the number of edges in shortest binary paths (i.e. shortest
%   paths computed on the binarized connectivity matrix), because shortest 
%   weighted paths have the minimal weighted distance, but not necessarily 
%   the minimal number of edges.
%       Lengths between disconnected nodes are set to Inf.
%       Lengths on the main diagonal are set to 0.
%
%   Algorithm: Dijkstra's algorithm.
%
%
%   Mika Rubinov, UNSW/U Cambridge, 2007-2012.
%   Rick Betzel and Andrea Avena, IU, 2012

%Modification history
%2007: original (MR)
%2009-08-04: min() function vectorized (MR)
%2012: added number of edges in shortest path as additional output (RB/AA)
%2013: variable names changed for consistency with other functions (MR)


%%
%CHARPATH       Characteristic path length, global efficiency and related statistics
%
%   lambda                                  = charpath(D);
%   lambda                                  = charpath(D);
%   [lambda,efficiency]                     = charpath(D);
%   [lambda,efficiency,ecc,radius,diameter] = charpath(D,diagonal_dist,infinite_dist);
%
%   The network characteristic path length is the average shortest path
%   length between all pairs of nodes in the network. The global efficiency
%   is the average inverse shortest path length in the network. The nodal
%   eccentricity is the maximal path length between a node and any other
%   node in the network. The radius is the minimal eccentricity, and the
%   diameter is the maximal eccentricity.
%
%   Input:      D,              distance matrix
%               diagonal_dist   optional argument
%                               include distances on the main diagonal
%                                   (default: diagonal_dist=0)
%               infinite_dist   optional argument
%                               include infinite distances in calculation
%                                   (default: infinite_dist=1)
%
%   Outputs:    lambda,         network characteristic path length
%               efficiency,     network global efficiency
%               ecc,            nodal eccentricity
%               radius,         network radius
%               diameter,       network diameter
%
%   Notes:
%       The input distance matrix may be obtained with any of the distance
%   functions, e.g. distance_bin, distance_wei.
%       Characteristic path length is defined here as the mean shortest
%   path length between all pairs of nodes, for consistency with common
%   usage. Note that characteristic path length is also defined as the
%   median of the mean shortest path length from each node to all other
%   nodes.
%       Infinitely long paths (i.e. paths between disconnected nodes) are
%   included in computations by default. This behavior may be modified with
%   via the infinite_dist argument.
%
%
%   Olaf Sporns, Indiana University, 2002/2007/2008
%   Mika Rubinov, U Cambridge, 2010/2015

%   Modification history
%   2002: original (OS)
%   2010: incorporation of global efficiency (MR)
%   2015: exclusion of diagonal weights by default (MR)
%   2016: inclusion of infinite distances by default (MR)

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

%% Set negative values to 0
% set negatives to 0

PatientTable3d_positive=PatientTable3d_nan_to_0;

PatientTable3d_positive(PatientTable3d_positive < 0) = 0;

%% Making sure everything is correct
%The matrix with r to z tranforms becomes your master_matrix

master_connectome(:,:,:)=PatientTable3d_positive;% Change master_connectome to connectome matrix

NP=size(master_connectome,3); % Number of participants.
AtlasSize=size(master_connectome,1); % Ensure this value is connectome size. Depends on atlas

%% Convert weighted matrix to length matrix
% Weighted matrix must be converted to length matrix. Length is then fed
% in to get distancce matrix script. 

%obtains 3d matrix for length matrices

for i =1:NP
    L(:,:,i)=master_connectome(:,:,i);
    E=find(L); %Finds every position for value in L does not equal 0
    L(E)=1./L(E);           % invert weights. L(E) is all values at position of E 1 divided by all values in L. 
    length_matrix(:,:,i)=L(:,:,i);
end


%% Convert L zeros to inf and then set diagonal to 0
%Set all 0s to inf to indicate no connection between nodes. 
length_matrix(length_matrix==0)=Inf;

%Then set diagonal to 0
for n = 1:size(length_matrix,1)
    length_matrix(n,n,:) = 0;
end

%% Convert length matrix to distance matrix
%CURRENTLY WORKING

for i = 1:NP
    %loop through 87 length matrices
    L_per_subj=length_matrix(:,:,i);
    n=length(L_per_subj); %length should be length of matrix
    D=inf(n); %make matrix filled with inf, diagonal set to 0
    D(1:n+1:end)=0;
    B=zeros(n);     %matrix of 0s               %number of edges matrix
    for u=1:n %row of 1:length of matrix. for ex, 1-400
        S=true(1,n);        %row filled with 1's, logical argument       %distance permanence (true is temporary)
        L1=L_per_subj;
        V=u; %Number of nodes, e.g., 400
        while 1
            S(V)=0; %Fill with 0s instead of 1s    %distance u->V is now permanent
            L1(:,V)=0; %Fill with 0s                         %no in-edges as already shortest
            for v=V % indeces per node (e.g., 400)
                T=find(L1(v,:));                %neighbours of shortest nodes
                [d,wi]=min([D(u,T);D(u,v)+L1(v,T)]); %find smallest value
                D(u,T)=d;                       %smallest of old/new path lengths
                ind=T(wi==2);                   %indices of lengthened paths
                B(u,ind)=B(u,v)+1;              %increment no. of edges in lengthened paths
            end

            minD=min(D(u,S));
            if isempty(minD)||isinf(minD)       %isempty: all nodes reached;
                break,                          %isinf: some nodes cannot be reached
            end;

            V=find(D(u,:)==minD);
            
            %Save the short path length (B) and distance matrixes (D) for
            %each participant!!
            distance_matrix(:,:,i)=D;
            b_short_edges(:,:,i)=B;

        end
    end
end

%% Find charpath 

%               lambda,         network characteristic path length
%               efficiency,     network global efficiency
%               ecc,            nodal eccentricity

%optional, I'm not using these but keeping for future reference. 
%               radius,         network radius
%               diameter,       network diameter

n = size(D,1);

for i=1:NP 
    % if any(any(isnan(distance_matrix(:,:,i))))
        % error('The distance matrix must not contain NaN values');
    % end
    
    for n = 1:size(distance_matrix,1)
    distance_matrix(n,n,:) = NaN; % set diagonal distance to NaN
    end          
    
    distance_matrix_nan=distance_matrix;

    distance_matrix_nan(isinf(distance_matrix_nan))  = NaN;             % ignore infinite path lengths
    
    distance_matrix_path=distance_matrix_nan(:,:,i);

    Dv = distance_matrix_path(~isnan(distance_matrix_path));                  % get non-NaN indices of D
    
    % Mean of entries of D(G)
    lambda(i,:)     = mean(Dv);
    
    % Efficiency: mean of inverse entries of D(G)
    efficiency(i,:) = mean(1./Dv);
    
    % Eccentricity for each vertex
    ecc(i,:)       = nanmax(distance_matrix_path,[],2);

    %Eccentricity mean (not included in original charpath script)
    ecc_avg(i,:)       = nanmean(ecc(i,:));
    
    % Radius of graph
    radius(i,:)     = min(ecc(i,:));
    
    % Diameter of graph
    diameter(i,:)   = max(ecc(i,:));
end

