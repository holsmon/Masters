function [S, W, B, S_Mod, W_Mod, B_Mod] = w_b_s_calc_with_nan(M, Ci) 
% DESCRIPTION:
%   The degree to which edges are more dense within communites and more 
%   sparse between communities, which quantifies the segregation of a 
%   weighted network (Chan et al. 2014). 
%
% USAGE: 
%   S = segregation(M,Ci)
%   [S, W, B] = segregation(M,Ci)
%
% Inputs:   M,      weighted symmetrical matrix (n x n matrix)
%           Ci,     community affiliation vector (n x 1 vector)
%
% Outputs:  S,      System segregation calcualted with W & B                
%           W,      mean edge weight between nodes within the same
%                   community     
%           B,      mean edge weight obetween nodes from different
%                   community  
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reference: Chan et al. (2014) PNAS E4997
%   2014-2018
%   Micaela Chan, UTD
%
%   Modification History:
%   2014: original
%   Oct 2018: commented script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Notes about diagonal (not meaningful since this represents a nodes connection to itself)

%This script is NOT including the diagonal, so the diagonal cell count is not 
%being counted in the denominator. 
%this is because the script only takes upper right triangle of the matrix,
%which ignores the diagonal. 
%This is an important distinction if you're trying to
%replicate the values this script is producing - HM 8/4/2021

%% Notes about missing nodes 
%- HM 5/10/2022

%some participants in VisAtten.03 had missing nodes due to slice cutoff.
%Therefore, nodes without data were set to NaN. I used the MATLAB function
%nanmean to account for this. 

%% Begin analysis

if length(Ci)~=length(M)
        error('Length of label does not match with matrix dimension.');
end

nCi = unique(Ci);
Wv = [];
Bv = [];
W_Mod = [];
B_Mod = [];
S_Mod = [];

%indexing what is within a module and what is between for 192 nodes. Gives
%1x192 matrix where each node is labelled 0 or 1 to indicate if its within
%or between. Does this for every single module.  
for i = 1:length(nCi) % loop through communities
   Wi = Ci == nCi(i); % find index for within communitiy edges
   Bi = Ci ~= nCi(i); % find index for between communitiy edges
   
   %extracts all within connectivity values for each module. 
   Wv_temp = M(Wi,Wi); % extract within communitiy edges
   Bv_temp = M(Wi,Bi); % extract between communitiy edges
      
   Wv_temp = Wv_temp(logical(triu(ones(sum(Wi)),1)))';
   %Wv is all within connections and Bv is all between communities. 
   Wv = [Wv, Wv_temp];  
   Bv = [Bv, Bv_temp(:)'];

   %Used nanmean so MATLAB wouldn't freak out at my nans in the
   %master_connectome - HM 5/10/2022
   % Get Subnetworks - MR
   %means Wv temp because that has connections per module. 
   W_Mod(1,i) = nanmean(Wv_temp);
   B_Mod(1,i) = nanmean(Bv_temp,'all');
end

%When calculating the mean, it's important to account for the fact that
%each module has a different number of nodes. Therefore, your mean should
%be weighted by the 7 different modules and how many nodes reside in each
%of the 7 modules. - HM 8/4/2021

% Global
W = nanmean(Wv); % mean within community edges 
B = nanmean(Bv); % mean between community edges
S = (W-B)/W; % system segregation

% Subnetwork System Segregation
S_Mod = ((W_Mod - B_Mod)./W_Mod).';


