function [vi_expr,vi_gen,vi_sum,vi_prod] = GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx,input_idx,tree_method,K,ntrees)
%Computation of tree-based weights for putative edges directed towards a 
%specified target gene.
%
%[vi_expr,vi_gen,vi_sum,vi_prod] = 
%GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx) learns a tree
%model from expr_matrix and gen_matrix and assigns a weight to each edge
%directed from a putative regulator to the target gene. 
%   * expr_matrix is a matrix containing expression values. Each line
%   corresponds to an experiment and each column corresponds to a gene. 
%   * gen_matrix is a matrix containing genotype data (only 0 and 1 values). 
%   Each line corresponds to an experiment and each column corresponds to 
%   a gene. Furthermore each line i (resp. column j) has to correspond to 
%   line i (resp. column j) of expr_matrix.
%   * output_idx is the (column) index of the target gene in expr_matrix
%   (and gen_matrix).
%vi_expr and vi_gen are vectors of length p, where p is the number of
%columns in expr_matrix (and gen_matrix). vi_expr(i) is the importance of
%the expression of gene i in the prediction of the expression of target
%gene. vi_gen(i) is the importance of the genotype of gene i in the
%prediction of the expression of target gene.
%vi_sum(i) = vi_expr(i) + vi_gen(i)
%vi_prod(i) = vi_expr(i) * vi_gen(i)
%
%[vi_expr,vi_gen,vi_sum,vi_prod] = 
%GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx,input_idx) only
%uses as candidate regulators the genes whose index (as ordered in
%expr_matrix and gen_matrix) is in input_idx. input_idx is a vector of
%length <= p. vi_expr(i), vi_gen(i), vi_sum(i) and vi_prod(i) such that i
%is not in input_idx are set to zero. The default vector contains the
%indices of all the genes in expr_matrix and gen_matrix.
%
%[vi_expr,vi_gen,vi_sum,vi_prod] = 
%GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx,input_idx,tree_me
%thod) specifies which tree precedure is used. Available methods are:
%   * 'RF' - Random Forests (Default method)
%   * 'ET' - Extra Trees
%
%[vi_expr,vi_gen,vi_sum,vi_prod] = 
%GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx,input_idx,tree_me
%thod,K) specifies the number K of randomly selected attributes at each
%node of one tree. Possible values of K:
%   * 'sqrt' - K = square root of the number of candidate regulators
%   (Default value)
%   * 'all' - K = number of candidate regulators
%   * any numerical value
%
%[vi_expr,vi_gen,vi_sum,vi_prod] =
%GENIE3_SG_joint_single(expr_matrix,gen_matrix,output_idx,input_idx,tree_me
%thod,K,ntrees) specifies the number of trees grown in the ensemble.
%Default value: 1000.


%% Check input arguments
error(nargchk(3,7,nargin));

if size(expr_matrix,1) ~= size(gen_matrix,1) || size(expr_matrix,2) ~= size(gen_matrix,2)
    error('expr_matrix and gen_matrix must be of same size.')
end

if length(output_idx) ~= 1
    error('Input argument output_idx must be one integer.')
end

if ~ismember(output_idx,1:size(expr_matrix,2))
    error('Input argument output_idx must be an integer between 1 and p, where p is the number of genes in expr_matrix and gen_matrix.')
end
   
if nargin > 3 && sum(ismember(input_idx,1:size(expr_matrix,2))) ~= length(input_idx)
    error('Input argument input_idx must be a vector containing integers between 1 and p, where p is the number of genes in expr_matrix.')
end

if nargin > 4 && sum(strcmp(tree_method,{'RF' 'ET'})) == 0
    error('Input argument tree_method must be ''RF'' or ''ET''.')
end

if nargin > 5 && isa(K,'char') && sum(strcmp(K,{'sqrt' 'all'})) == 0
    error('Input argument K must be ''sqrt'', ''all'' or a numerical value.')
end

if nargin > 6 && ~isa(ntrees,'numeric')
    error('Input argument ntrees must be an integer.')
end

%% Data must be in single precision when used to build a tree model
expr_matrix = single(expr_matrix);
gen_matrix = single(gen_matrix);
nsamples = size(expr_matrix,1); % number of experiments
ngenes = size(expr_matrix,2); % number of genes

%% Output vector
output = expr_matrix(:,output_idx);
% Normalize output data
output = output/std(output);


%% Indices of input genes
if nargin >= 4
    input_idx = unique(input_idx);
else
    % Default: all the genes are candidate regulators
    input_idx = 1:ngenes;
end

% Remove target gene from candidate regulators (for the expression data)
input_idx_expr = setdiff(input_idx,output_idx);

% Add target gene to candidate regulators (for the genotype data)
if ~ismember(output_idx,input_idx)
    input_idx_gen = [input_idx output_idx];
else
    input_idx_gen = input_idx;
end

ninputs = length(input_idx_expr) + length(input_idx_gen);


%% Tree parameters

% Default parameters: Random Forests, K=sqrt(number of input genes),
% 1000 trees in the ensemble
if nargin < 5 || (nargin >= 5 && strcmp(tree_method,'RF'))
    ok3ensparam = init_rf();
    if nargin >= 6
        if strcmp(K,'all')
            ok3ensparam=init_rf(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_rf(round(K));
        end
    end
elseif nargin >= 5 && strcmp(tree_method,'ET')
    ok3ensparam = init_extra_trees();
    if nargin >= 6
        if strcmp(K,'all')
            ok3ensparam=init_extra_trees(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_extra_trees(round(K));
        end
    end
end

% Number of trees in the ensemble
if nargin < 7
    ok3ensparam.nbterms = 1000;
else
    ok3ensparam.nbterms = ntrees;
end
    

%% Learning of tree model
[tree,varimp]=rtenslearn_c([expr_matrix(:,input_idx_expr) gen_matrix(:,input_idx_gen)],output,int32(1:nsamples),[],ok3ensparam,0);
vi_expr = zeros(1,ngenes);
vi_gen = zeros(1,ngenes);
varimp = varimp';
vi_expr(input_idx_expr) = varimp(1:length(input_idx_expr)) / nsamples;
vi_gen(input_idx_gen) = varimp(length(input_idx_expr)+1:ninputs) / nsamples;
vi_expr(output_idx) = 0;
vi_gen(output_idx) = 0;
vi_sum = vi_expr + vi_gen;
vi_prod = vi_expr .* vi_gen;
