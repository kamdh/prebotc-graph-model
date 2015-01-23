import numpy as np

def fit_MRF_pseudolikelihood(adj_exc,adj_inh,y):
    '''
    Fit a Markov random field using maximum pseudolikelihood estimation,
    also known as logistic regression. The conditional probabilities
    follow

    y_i ~ Logistic(B[0] + B[1] A1_{ij} y_j + A1[2] X_{ij} (1-y_j)
                   + B[3] A2_{ij} y_j + B[4] A3_{ij} (1-y_j) ),

    where A1 = adj_exc and A2 = adj_inh and each term is summed over
    j.

    Params
    ======
      adj_exc: excitatory adjacency matrix
      adj_inh: inhibitory adjacency matrix
      y: site variables, 0 or 1

    Returns
    =======
      B: logistic regression coefficients
    '''
    from sklearn.linear_model import LogisticRegression
    from sklearn import cross_validation
    if len(np.unique(y)) < 2:
        B=np.array([np.nan,np.nan,np.nan,np.nan,np.nan])
    else:
        N=y.shape[0]
        ytile=np.tile(y,(N,1)).T
        X1=np.array(np.sum(np.multiply(adj_exc,ytile),0)).flatten()
        X2=np.array(np.sum(np.multiply(adj_exc,1-ytile),0)).flatten()
        X3=np.array(np.sum(np.multiply(adj_inh,ytile),0)).flatten()
        X4=np.array(np.sum(np.multiply(adj_inh,1-ytile),0)).flatten()
        model=LogisticRegression(penalty='l2')
        X=np.column_stack((X1,X2,X3,X4))
        model.fit(X,y)
        B=model.raw_coef_.flatten()
    return B
    
def predict_MRF(B, adj_exc, adj_inh, burn_in=4e3, steps=1e4,
                skip_multiple=3):
    '''
    Perform prediction with an MRF (Markov random field). Uses Gibbs sampling
    to sample from the distribution P(y) =1/Z exp( -H(y) ).

    The Hamiltonian is:

    H = \sum_{ij} y_i (B[0] + B[1] A1_{ji} y_j + A1[2] X_{ji} (1-y_j)
                       + B[3] A2_{ji} y_j + B[4] A3_{ji} (1-y_j))

    Params
    ======
      B: coefficients of the MRF
      adj_exc: excitatory adjacency matrix
      adj_inh: inhibitory adjacency matrix
      burn_in: number of burn-in steps to take (default: 4000)
      steps: total number of Gibbs steps to take (default: 10000)
      skip_multiple: skips skip_multiple * num_neuron steps between samples

    Returns
    =======
      ypostmean: posterior mean of state
    '''
    import numpy.random
    def gibbs_proba(y,B,adj_exc,adj_inh):
        term0=B[0]*np.dot(adj_exc.T,y)
        term1=B[1]*np.dot(adj_exc.T,(1-y))
        term2=B[2]*np.dot(adj_inh.T,y)
        term3=B[3]*np.dot(adj_inh.T,(1-y))
        e=B[4]+term0+term1+term2+term3
        return np.exp(e)/(np.exp(e)+1.0)
    
    N=adj_exc.shape[0]
    steps=int(steps)
    # run a Gibbs sampler
    y=np.random.rand(N,1)
    samples=np.zeros((N,steps))
    # zero diagonals (should be 0 already)
    np.fill_diagonal(adj_exc,0) 
    np.fill_diagonal(adj_inh,0)
    for ii in range(steps):
        yt=y
        # Gibbs update
        proba=gibbs_proba(y,B,adj_exc,adj_inh)
        yt=np.array(np.random.rand(N,1) < proba,dtype=np.float)
        y=yt
        samples[:,ii]=y.flatten()
    # compute mle
    use_indices=np.arange(burn_in,steps,skip_multiple*N, dtype=int)
    final_samples=samples[:,use_indices]
    ypostmean=np.mean(final_samples,axis=1)
    return ypostmean

def fit_logistic_graph_features():
    pass

def get_node_features(adj_exc,adj_inh,normalize_centrality=True):
    '''
    Get node-based features to train logistic classifier
    
    Params
    ======
      adj_exc: excitatory adjacency matrix
      adj_inh: inhibitory adjacency matrix
      normalize_centrality: normalize relevant measures? (default: True)

    Returns
    =======
      X: numneuron x numfeatures array to be used with logistic regression
      X_labels
    '''
    import networkx as nx
    G_exc=nx.DiGraph(adj_exc)
    G_inh=nx.DiGraph(adj_inh)
    
    def dict_to_array(d):
        return np.array([d[i] for i in sorted(d)])

    def features(G,normalize_centrality):
        '''
        Returns the features we are interested in within a dict
        '''
        load_centrality=nx.load_centrality(G,normalized=normalize_centrality)
        betweenness_centrality=nx.betweenness_centrality(G,normalized=normalize_centrality)
        eigenvector_centrality=nx.eigenvector_centrality_numpy(G,normalized=normalize_centrality)
        closeness_centrality=nx.closeness_centrality(G,normalized=normalize_centrality)
        in_degree=G.in_degree()
        out_degree=G.out_degree()
        core_number=nx.core_number(G)
        clustering=nx.clustering(G)
        d={}
        d['in_degree']=in_degree
        d['out_degree']=out_degree
        d['load_centrality']=load_centrality
        d['betweennes_centrality']=betweennes_centrality
        d['eigenvector_centrality']=eigenvector_centrality
        d['closeness_centrality']=closeness_centrality
        d['core_number']=core_number
        return d

    # grab the features
    d_exc=features(G_exc)
    d_inh=features(G_inh)
    # setup some structures
    num_features=len(d_exc)+len(d_inh)
    num_nodes=G.number_of_nodes()
    X=np.zeros((num_nodes,num_features),dtype=np.float)
    X_labels=[]
    # fill in X and Xlabels
    feature_index=0
    for gclass in ('exc','inh'):
        if gclass == 'exc':
            d=d_exc
        else:
            d=d_inh
        for feature in sorted(d):
            X_labels.append(feature+"_"+gclass)
            X[:,feature_index]=dict_to_array(d[feature])
            feature_index+=1
    return X, X_labels
