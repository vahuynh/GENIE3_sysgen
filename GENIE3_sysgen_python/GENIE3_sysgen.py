from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.tree import BaseDecisionTree
from numpy import *
import time
from operator import itemgetter
from multiprocessing import Pool

def compute_feature_importances(estimator):
    
    """Computes variable importances from a trained tree-based model.
    """
    
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = array(importances)
        return sum(importances,axis=0) / len(estimator)


def get_link_list(VIM,gene_names=None,regulators='all',maxcount='all',file_name=None):
    
    """Gets the ranked list of (directed) regulatory links.
    
    Parameters
    ----------
    
    VIM: numpy array
        Array as returned by the function GENIE3(), in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. 
        
    gene_names: list of strings, optional
        List of length p, where p is the number of rows/columns in VIM, containing the names of the genes. The i-th item of gene_names must correspond to the i-th row/column of VIM. When the gene names are not provided, the i-th gene is named Gi.
        default: None
        
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names), and the returned list contains only edges directed from the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    maxcount: 'all' or positive integer, optional
        Writes only the first maxcount regulatory links of the ranked list. When maxcount is set to 'all', all the regulatory links are written.
        default: 'all'
        
    file_name: string, optional
        Writes the ranked list of regulatory links to the file file_name.
        default: None
        
        
    
    Returns
    -------
    
    The list of regulatory links, ordered according to the edge score. Auto-regulations do not appear in the list. Regulatory links with a score equal to zero are randomly permuted. In the ranked list of edges, each line has format:
        
        regulator   target gene     score of edge
    """
    
    # Check input arguments      
    if not isinstance(VIM,ndarray):
        raise ValueError('VIM must be a square array')
    elif VIM.shape[0] != VIM.shape[1]:
        raise ValueError('VIM must be a square array')
        
    ngenes = VIM.shape[0]
        
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')
        
    if maxcount is not 'all' and not isinstance(maxcount,int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')
        
    if file_name is not None and not isinstance(file_name,str):
        raise ValueError('input argument file_name must be a string')
    
    

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
               
    nTFs = len(input_idx)
    
    # Get the non-ranked list of regulatory links
    vInter = [(i,j,score) for (i,j),score in ndenumerate(VIM) if i in input_idx and i!=j]
    
    # Rank the list according to the weights of the edges        
    vInter_sort = sorted(vInter,key=itemgetter(2),reverse=True)
    nInter = len(vInter_sort)
    
    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx,target_idx,score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1
            
    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = random.permutation(items_perm)
        vInter_sort[i:] = items_perm
        
    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount,int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount
        
    if file_name:
    
        outfile = open(file_name,'w')
    
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx+1,target_idx+1,score))
            
        
        outfile.close()
        
    else:
        
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print '%s\t%s\t%.6f' % (gene_names[TF_idx],gene_names[target_idx],score)
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print 'G%d\tG%d\t%.6f' % (TF_idx+1,target_idx+1,score)
                
                
    
    
def GENIE3(expr_data,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,nthreads=1):
    
    '''Computation of tree-based scores for all putative regulatory links.
    
    Parameters
    ----------
    
    expr_data: numpy array
        Array containing gene expression values. Each row corresponds to a condition and each column corresponds to a gene.
        
    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data, containing the names of the genes. The i-th item of gene_names must correspond to the i-th column of expr_data.
        default: None
        
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
        
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
         
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000
    
    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
        
        
    Returns
    -------
    
    An array in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. All diagonal elements are set to zero (auto-regulations are not considered). When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.
    
    '''
    
    time_start = time.time()
    
    # Check input arguments
    if not isinstance(expr_data,ndarray):
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
        
    ngenes = expr_data.shape[1]
    
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expr_data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('the genes must contain at least one candidate regulator')        
        
    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
        
    if K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
        
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    
    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
        
    if not isinstance(nthreads,int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
        
    
    print 'Tree method: ' + str(tree_method)
    print 'K: ' + str(K)
    print 'Number of trees: ' + str(ntrees)
    print '\n'
    
    
    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]

    
    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    VIM = zeros((ngenes,ngenes))

    if nthreads > 1:
        print 'running jobs on %d threads' % nthreads

        input_data = list()
        for i in range(ngenes):
            input_data.append( [expr_data,i,input_idx,tree_method,K,ntrees] )

        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_single, input_data)
    
        for (i,vi) in alloutput:
            VIM[i,:] = vi

    else:
        print 'running single threaded jobs'
        for i in range(ngenes):
            print 'Gene %d/%d...' % (i+1,ngenes)
            
            VIM[i,:] = GENIE3_single(expr_data,i,input_idx,tree_method,K,ntrees)
   
    VIM = transpose(VIM)
   
    time_end = time.time()
    print "Elapsed time: %.2f seconds" % (time_end - time_start)

    return VIM
    
    
    
def wr_GENIE3_single(args):
    return([args[1], GENIE3_single(args[0], args[1], args[2], args[3], args[4], args[5])])
    


def GENIE3_single(expr_data,output_idx,input_idx,tree_method,K,ntrees):
    
    ngenes = expr_data.shape[1]
    
    # Expression of target gene
    output = expr_data[:,output_idx]
    
    # Normalize output data
    output = output / std(output)
    
    # Remove target gene from candidate regulators
    input_idx = input_idx[:]
    if output_idx in input_idx:
        input_idx.remove(output_idx)

    expr_data_input = expr_data[:,input_idx]

    # Parameter K of the tree-based method
    if (K == 'all') or (isinstance(K,int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K
    
    if tree_method == 'RF':
        treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=max_features)
    elif tree_method == 'ET':
        treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=max_features)

    # Learn ensemble of trees
    treeEstimator.fit(expr_data_input,output)
    
    # Compute importance scores
    feature_importances = compute_feature_importances(treeEstimator)
    vi = zeros(ngenes)
    vi[input_idx] = feature_importances
    
    return vi
        
    
    
    
    
    
def GENIE3_mark(expr_data,gen_data,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,nthreads=1):
    
    '''Computation of tree-based scores for all putative regulatory links directed from genotypes to gene expressions.
    
    Parameters
    ----------
    
    expr_data: numpy array
        Array containing gene expression values. 
        Each row corresponds to a condition and each column corresponds to a gene.
        
    gen_data: numpy array
        Array containing genotype values (only 0/1 values). 
        Each row corresponds to a condition and each column corresponds to a gene. Each line i (resp. column j) must correspond to 
        line i (resp. column j) of expr_data.
        
    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data and gen_data, containing the names of the genes. 
        The i-th item of gene_names must correspond to the i-th column of expr_data and gen_data.
        default: None
        
    regulators: list of string, optional
        List containing the names of the candidate regulators. 
        When a list of regulators is provided, the names of all the genes must be provided (in gene_names). 
        When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
        
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
         
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000
    
    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
        
        
    Returns
    -------
    
    An array in which the element (i,j) is the score of the edge directed from the genotype of the i-th gene to the expression of the j-th gene. 
    All diagonal elements are set to zero (auto-regulations are not considered). 
    When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.
    
    '''
    
    time_start = time.time()

    
    # Check input arguments
    if not isinstance(expr_data,ndarray):
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
        
    if not isinstance(gen_data,ndarray):
        raise ValueError('gen_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
    
    nsamples = expr_data.shape[0]    
    ngenes = expr_data.shape[1]
    
    if gen_data.shape[0]!= nsamples or gen_data.shape[1] != ngenes:
        raise ValueError('expr_data and gen_data must have the same numbers of rows and columns')
    
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')        
        
    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
        
    if K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
        
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    
    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
        
    if not isinstance(nthreads,int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
        
        
    print 'Tree method: ' + str(tree_method)
    print 'K: ' + str(K)
    print 'Number of trees: ' + str(ntrees)
    print '\n'
    
    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]

    
    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    VIM = zeros((ngenes,ngenes))
    
    if nthreads > 1:
        print 'running jobs on %d threads' % nthreads
        
        input_data = list()
        for i in range(ngenes):
            input_data.append( [expr_data,gen_data,i,input_idx,tree_method,K,ntrees] )

        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_mark_single, input_data)
    
        for (i,vi) in alloutput:
            VIM[i,:] = vi
            
    else:
        print 'running single threaded jobs'
        for i in range(ngenes):
            print 'Gene %d/%d...' % (i+1,ngenes)
            
            VIM[i,:] = GENIE3_mark_single(expr_data,gen_data,i,input_idx,tree_method,K,ntrees)
        
    VIM = transpose(VIM)
    
    time_end = time.time()
    print "Elapsed time: %.2f seconds" % (time_end - time_start)

    return VIM
    
    
def wr_GENIE3_mark_single(args):
    return([args[2], GENIE3_mark_single(args[0], args[1], args[2], args[3], args[4], args[5], args[6])])

    
    
def GENIE3_mark_single(expr_data,gen_data,output_idx,input_idx,tree_method,K,ntrees):

    ngenes = gen_data.shape[1]
    
    # Expression of target gene
    output = expr_data[:,output_idx]
    
    # Normalize output data
    output = output / std(output)

    # Add target gene to candidate regulators
    input_idx = input_idx[:]
    if output_idx not in input_idx:
        input_idx.append(output_idx)
        
    # Input data 
    matrix_input = gen_data[:,input_idx]

    # Parameter K of the tree-based method
    if (K == 'all') or (isinstance(K,int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K
    
    if tree_method == 'RF':
        treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=max_features)
    elif tree_method == 'ET':
        treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=max_features)

    # Learn ensemble of trees
    treeEstimator.fit(matrix_input,output)

    # Compute importance scores
    feature_importances = compute_feature_importances(treeEstimator)
    vi = zeros(ngenes)
    vi[input_idx] = feature_importances
    vi[output_idx] = 0

    return vi
    
    
    
    
    
def GENIE3_SG_sep(expr_data,gen_data,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,nthreads=1):
    
    '''Computation of tree-based scores for all putative regulatory links. For each gene, two separate tree-based models are learned,
    from the expression and genotype data respectively.
    
    Parameters
    ----------
    
    expr_data: numpy array
        Array containing gene expression values. 
        Each row corresponds to a condition and each column corresponds to a gene.
        
    gen_data: numpy array
        Array containing genotype values (only 0/1 values). 
        Each row corresponds to a condition and each column corresponds to a gene. Each line i (resp. column j) must correspond to 
        line i (resp. column j) of expr_data.
        
    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data and gen_data, containing the names of the genes. 
        The i-th item of gene_names must correspond to the i-th column of expr_data and gen_data.
        default: None
        
    regulators: list of string, optional
        List containing the names of the candidate regulators. 
        When a list of regulators is provided, the names of all the genes must be provided (in gene_names). 
        When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
        
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
         
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000
    
    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
        
        
    Returns
    -------
    
    (VIM_expr, VIM_gen, VIM_sum, VIM_prod): 4 arrays of size p x p, where p is the number of genes.
    
    VIM_expr[i,j] is the score of the expression of gene i in the prediction of the expression of gene j.
    
    VIM_gen[i,j] is the score of the genotype of gene i in the prediction of the expression of gene j.
    
    VIM_sum[i,j] = VIM_expr[i,j] + VIM_gen[i,j]
    
    VIM_prod[i,j] = VIM_expr[i,j] * VIM_gen[i,j]
     
    All diagonal elements are set to zero (auto-regulations are not considered). 
    When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.
    
    '''

    
    time_start = time.time()
    
    # Check input arguments
    if not isinstance(expr_data,ndarray):
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
        
    if not isinstance(gen_data,ndarray):
        raise ValueError('gen_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
    
    nsamples = expr_data.shape[0]    
    ngenes = expr_data.shape[1]
    
    if gen_data.shape[0]!= nsamples or gen_data.shape[1] != ngenes:
        raise ValueError('expr_data and gen_data must have the same numbers of rows and columns')
    
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')        
        
    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
        
    if K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
        
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    
    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
        
    if not isinstance(nthreads,int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
        
    
    print '\nLearning tree models from expression data...'
    VIM_expr = GENIE3(expr_data,gene_names=gene_names,regulators=regulators,tree_method=tree_method,K=K,ntrees=ntrees,nthreads=nthreads)
    
    print '\nLearning tree models from genotype data...'
    VIM_gen = GENIE3_mark(expr_data,gen_data,gene_names=gene_names,regulators=regulators,tree_method=tree_method,K=K,ntrees=ntrees,nthreads=nthreads)
    
    VIM_sum = VIM_expr + VIM_gen
    VIM_prod = VIM_expr * VIM_gen
    
    return VIM_expr,VIM_gen,VIM_sum,VIM_prod
    
    
    
    
def GENIE3_SG_joint(expr_data,gen_data,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,nthreads=1):
    
    '''Computation of tree-based scores for all putative regulatory links. For each gene, one tree-based model is learned,
    jointly from the expression and genotype data.
    
    Parameters
    ----------
    
    expr_data: numpy array
        Array containing gene expression values. 
        Each row corresponds to a condition and each column corresponds to a gene.
        
    gen_data: numpy array
        Array containing genotype values (only 0/1 values). 
        Each row corresponds to a condition and each column corresponds to a gene. Each line i (resp. column j) must correspond to 
        line i (resp. column j) of expr_data.
        
    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data and gen_data, containing the names of the genes. 
        The i-th item of gene_names must correspond to the i-th column of expr_data and gen_data.
        default: None
        
    regulators: list of string, optional
        List containing the names of the candidate regulators. 
        When a list of regulators is provided, the names of all the genes must be provided (in gene_names). 
        When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
        
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
         
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000
    
    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
        
        
    Returns
    -------
    
    (VIM_expr, VIM_gen, VIM_sum, VIM_prod): 4 arrays of size p x p, where p is the number of genes.
    
    VIM_expr[i,j] is the score of the expression of gene i in the prediction of the expression of gene j.
    
    VIM_gen[i,j] is the score of the genotype of gene i in the prediction of the expression of gene j.
    
    VIM_sum[i,j] = VIM_expr[i,j] + VIM_gen[i,j]
    
    VIM_prod[i,j] = VIM_expr[i,j] * VIM_gen[i,j]
     
    All diagonal elements are set to zero (auto-regulations are not considered). 
    When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.
    
    '''
    
    time_start = time.time()
    
    # Check input arguments
    if not isinstance(expr_data,ndarray):
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
        
    if not isinstance(gen_data,ndarray):
        raise ValueError('gen_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')
    
    nsamples = expr_data.shape[0]    
    ngenes = expr_data.shape[1]
    
    if gen_data.shape[0]!= nsamples or gen_data.shape[1] != ngenes:
        raise ValueError('expr_data and gen_data must have the same numbers of rows and columns')
    
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')        
        
    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
        
    if K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
        
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    
    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
        
    if not isinstance(nthreads,int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
        
    
    print 'Tree method: ' + str(tree_method)
    print 'K: ' + str(K)
    print 'Number of trees: ' + str(ntrees)
    print '\n'
        
        
    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
        
    
    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    VIM_expr = zeros((ngenes,ngenes))
    VIM_gen = zeros((ngenes,ngenes))
    VIM_sum = zeros((ngenes,ngenes))
    VIM_prod = zeros((ngenes,ngenes))
    
    
    if nthreads > 1:
        print 'running jobs on %d threads' % nthreads

        input_data = list()
        for i in range(ngenes):
            input_data.append( [expr_data,gen_data,i,input_idx,tree_method,K,ntrees] )

        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_SG_joint_single, input_data)
    
        for out in alloutput:
            i = out[0]
            
            (vi_expr,vi_gen,vi_sum,vi_prod) = out[1]
            VIM_expr[i,:] = vi_expr
            VIM_gen[i,:] = vi_gen
            VIM_sum[i,:] = vi_sum
            VIM_prod[i,:] = vi_prod
            
    else:
        print 'running single threaded jobs'
        for i in range(ngenes):
            print 'Gene %d/%d...' % (i+1,ngenes)
            
            (vi_expr,vi_gen,vi_sum,vi_prod) = GENIE3_SG_joint_single(expr_data,gen_data,i,input_idx,tree_method,K,ntrees)
            VIM_expr[i,:] = vi_expr
            VIM_gen[i,:] = vi_gen
            VIM_sum[i,:] = vi_sum
            VIM_prod[i,:] = vi_prod
    
        
    VIM_expr = transpose(VIM_expr)
    VIM_gen = transpose(VIM_gen)
    VIM_sum = transpose(VIM_sum)
    VIM_prod = transpose(VIM_prod)
    
    time_end = time.time()
    print "Elapsed time: %.2f seconds" % (time_end - time_start)

    return VIM_expr,VIM_gen,VIM_sum,VIM_prod



def wr_GENIE3_SG_joint_single(args):
    return([args[2], GENIE3_SG_joint_single(args[0], args[1], args[2], args[3], args[4], args[5], args[6])])
    
    
    
def GENIE3_SG_joint_single(expr_data,gen_data,output_idx,input_idx,tree_method,K,ntrees):

    nsamples = expr_data.shape[0]
    ngenes = expr_data.shape[1]

    # Expression of target gene
    output = expr_data[:,output_idx]
    
    # Normalize output data
    output = output / std(output)

    input_idx_expr = input_idx[:]
    input_idx_gen = input_idx[:]

    # Remove target gene from candidate regulators (for the expression data)
    if output_idx in input_idx_expr:
        input_idx_expr.remove(output_idx)
        
    # Add target gene to candidate regulators (for the genotype data)
    if output_idx not in input_idx_gen:
        input_idx_gen.append(output_idx)

    ninputs = len(input_idx_expr) + len(input_idx_gen)

    matrix_input = zeros((nsamples,ninputs))
    matrix_input[:,:len(input_idx_expr)] = expr_data[:,input_idx_expr]
    matrix_input[:,len(input_idx_expr):] = gen_data[:,input_idx_gen]


    # Parameter K of the tree-based method
    if (K == 'all') or (isinstance(K,int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K
    
    if tree_method == 'RF':
        treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=max_features)
    elif tree_method == 'ET':
        treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=max_features)
    
    
    # Learn ensemble of trees
    treeEstimator.fit(matrix_input,output)
    
    # Compute importance scores
    feature_importances = compute_feature_importances(treeEstimator)
    
    vi_expr = zeros(ngenes)
    vi_gen = zeros(ngenes)
    vi_expr[input_idx_expr] = feature_importances[:len(input_idx_expr)]
    vi_gen[input_idx_gen] = feature_importances[len(input_idx_expr):]
     
    vi_expr[output_idx] = 0
    vi_gen[output_idx] = 0
                
    vi_sum = vi_expr + vi_gen
    vi_prod = vi_expr * vi_gen


    return vi_expr,vi_gen,vi_sum,vi_prod
    

