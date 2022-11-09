
# good writeup
# https://www.multibugs.org/documentation/latest/spatial/SpatialDistributions.html#Appendix



# in REPL or VSCODE, this needs to be loaded first as "startup.jl" is skipped
project_directory = @__DIR__() #  same folder as the current file
push!(LOAD_PATH, project_directory)  # add the directory to the load path, so it can be found
include( "startup.jl" )


project_directory = "/home/jae/projects/dynamical_model/car"

 
using Turing 
using PDMats
using LinearAlgebra 
using StatsModels
using DataFrames
using SparseArrays
using Graphs


include( joinpath( project_directory, "car_functions.jl"  ))  


# ---------------
# load libs and options and prepare data for diffeq/turing model and set default parameters
include( joinpath( project_directory, "scottish_lip_cancer_environment.jl"  ))  # bootstrap different project environments depending on above choices


 



m = turing_car(D, W, X, log_offset, y)        # 204 sec

m = turing_car_prec(D, W, X, log_offset, y)   # 65 sec

m = turing_icar_direct_test( node1, node2, std(y) )  # Morris' "simple_iar" testing difference formulation .. using same run specs: results are similar with much better ess than stan   

m = turing_icar_direct_bym(X, log_offset, y, node1, node2) # 13 sec
   
# bym2 requires a "scaling factor"
m = turing_icar_direct_bym2(X, log_offset, y, node1, node2, scaling_factor_bym2(W)) # W is adjacency matrix , 18 sec; 30min for full run


# bym2 group model (multiple groups or disconnected groups): this is not finished 
scaling_factor = scaling_factor_bym2_groups(node1, node2, groups)
m = turing_icar_direct_bym2_groups(X, log_offset, y, node1, node2, scaling_factor, groups)


# use NUTS: see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz

# check timings and accuracy
n_samples, n_adapts, n_chains = 100, 100, 1
target_acceptance, max_depth, init_ϵ = 0.65, 7, 0.05
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)
o = sample(m, turing_sampler, n_samples) 


# larger runs (linux)
n_samples, n_adapts, n_chains = 9_000, 1_000, 4
target_acceptance, max_depth, init_ϵ = 0.65, 10, 0.1   # Morris uses 0.97 for target_acceptance, stan default is 0.95; such high acceptance rate does not work well -- divergent chains
turing_sampler = Turing.NUTS(n_adapts, target_acceptance; max_depth=max_depth, init_ϵ=init_ϵ)
o = sample( m, turing_sampler, MCMCThreads(), n_samples, n_chains  ) # to see progress
# o = sample(m, turing_sampler, n_samples) 

# if on windows and threads are still not working, use single processor mode:
# o = mapreduce(c -> sample(m, turing_sampler, n_samples), chainscat, 1:n_chains)

 

showall( summarize(o) )

 




int n;    // no. observations
int<lower=1> k; // no. of groups
int group_size[k]; // observational units per group
int group_idx[n]; // index of observations, ordered by group

int<lower=0> m; // no of components requiring additional intercepts
matrix[n, m] A; // dummy variables for any extra graph component intercepts
 
 
 
function icar_normal( 
    vector phi, 
    real spatial_scale,
    int[] node1, 
    int[] node2,
    int k, 
    int[] group_size, 
    int[] group_idx,
    int has_theta
    ) 
    
    real lp;
    int pos=1;

    # lp = -0.5 * dot_self( phi[node1] - phi[node2] );
    dphi = phi[node1] - phi[node2]
    lp = -0.5 * dot( dphi, dphi );

    if has_theta 
        for j in 1:k 
            # /* sum to zero constraint for each connected group; singletons zero out */
            # ig = segment(group_idx, pos, group_size[j])
            ig  = group_idx[ pos:(pos+group_size[j]-1) ]
            # lp += normal_lpdf( sum(phi[ ig ]) | 0, 0.001 * group_size[j] );
            Turing.@addlogprob! Normal( 0.0, 0.001 * group_size[j] )
            pos += group_size[j];
        end
    else 
        # does not have theta */
        for j in 1:k 
            if group_size[j] > 1 
                # /* same as above for non-singletons: sum to zero constraint */
                ig = group_idx[ pos:(pos+group_size[j]-1) ]
                # ig = segment(group_idx, pos, group_size[j])

                # lp += normal_lpdf(sum(phi[ig ]) | 0, 0.001 * group_size[j]);
                Turing.@addlogprob! Normal( 0.0, 0.001 * group_size[j] )
            else 
                # /* its a singleton: independent Gaussian prior on phi */
                ig = group_idx[ pos:(pos+group_size[j]-1) ]
                # ig = segment(group_idx, pos, group_size[j])
                # lp += normal_lpdf(phi[ ig  ] | 0, spatial_scale);
                Turing.@addlogprob! Normal( 0.0, spatial_scale )
            end
            pos += group_size[j];
        end
    end
    return lp;
end

 
# Combine local and global partial-pooling components into the convolved BYM2 term.
#  phi_tilde local (spatially autocorrelated) component
#  theta_tilde global component
#  spatial_scale scale parameter for the convolution term
#  n number of spatial units
#  k number of connected groups
#  group_size number of observational units in each group
#  group_idx index of observations in order of their group membership
#  inv_sqrt_scale_factor The scaling factor for the ICAR variance (see scale_c R function, using R-INLA);
#   transformed from 1/scale^2 --> scale. Or, a vector of ones.
#  rho proportion of convolution that is spatially autocorrelated
# # @return BYM2 convolution vector
 
# Scaling the ICAR prior

# To follow Reibler et al.’s adjustment to the scale of the model, you can use the INLA R package and the following R code:
 



library(INLA);
library(spdep);


# validateNb
# check that nbObject is symmetric, has connected components
#
validateNb = function(x) {
  if (is.symmetric.nb(x) && n.comp.nb(x)[[1]] < length(x)) return(TRUE);
  return(FALSE);
}


# isDisconnected
# check that nbObject has more than 1 component
#
isDisconnected = function(x) {
  return(n.comp.nb(x)[[1]] > 1);
}


# indexByComponent
#
# input: vector of component ids
# returns: vector of per-component consecutive node ids
#
indexByComponent = function(x) {
  y = x;
  comps = as.matrix(table(x));
  num_comps = nrow(comps);
  for (i in 1:nrow(comps)) {
    idx = 1;
    rel_idx = 1;
    while (idx <= length(x)) {
      if (x[idx] == i) {
        y[idx] = rel_idx;
        rel_idx = rel_idx + 1;
      }
      idx = idx + 1;
    }
  }
  return(y);
}


# nb2graph
#
# input: nb_object
# returns: dataframe containing num nodes, num edges,
#          and a list of graph edges from node1 to node2.
#
nb2graph = function(x) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (x[[i]][1] != 0) {
      n_links = n_links + length(x[[i]]);
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (x[[i]][1] > 0) {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i;
          node2[idx] = n2;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}


# nb2subgraph
# for a given subcomponent, return graph as lists of node1, node2 pairs
#
# inputs:
# x: nb object
# c_id: subcomponent id
# comp_ids: vector of subcomponent ids
# offsets: vector of subcomponent node numberings
# returns: list of node1, node2 ids
#
nb2subgraph = function(x, c_id, comp_ids, offsets) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        n_links = n_links + length(x[i]);
      }
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        for (j in 1:length(x[[i]])) {
          n2 = unlist(x[[i]][j]);    
          if (i < n2) {
            idx = idx + 1;
            node1[idx] = offsets[i];
            node2[idx] = offsets[n2];
          }
        }
      }
    }
  }
  return (list("node1"=node1,"node2"=node2));
}


# orderByComp
# given nbObject, reorder nb object so that all components are contiguous
# singletons moved to end of list
# returns list containing:
#   - new nbObject
#   - vector which maps old index to new index
#
orderByComponent = function(x) {
  if (!validateNb(x)) return(list());
  N = length(x);
  if (!isDisconnected(x)) return(list(x,seq(1:N)));

  rMap = rep(integer(0),N);
  comp_ids = n.comp.nb(x)[[2]];
  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  comp_sizes = as.vector(comps[,1]);
  idx = 1;
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] > 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  for (i in 1:nrow(comps)) {
    if (comp_sizes[i] == 1) {
      positions = which(comp_ids == i);
      for (j in 1:length(positions)) {
        rMap[idx] = as.integer(positions[j])
        idx = idx + 1;
      }
    }
  }
  new_ids = vector("character", length=N);
  for (i in 1:N) {
    idx_old = rMap[i];
    new_ids[i] = attributes(x)$region.id[idx_old];
  }

  # generate new nb list
  new_nb = structure(vector("list", length=N),class="nb");
  attr(new_nb, "region.id") = new_ids;  
  attr(new_nb, "type") = attributes(x)$type;
  attr(new_nb, "sym") = attributes(x)$sym;
  attr(new_nb, "region.id") = new_ids;  
  for (i in 1:N) {
    idx_old = rMap[i];
    old_nbs = x[[idx_old]];
    num_nbs = length(old_nbs);
    new_nb[[i]] = vector("integer", length=num_nbs);
    for (j in 1:num_nbs) {
      old_id = old_nbs[j];
      if (old_id == 0) {
        new_nb[[i]][j] = as.integer(0);
      } else {
        new_id = which(rMap == old_id);
        new_nb[[i]][j] = as.integer(new_id);
      }
    }
  }
  return(list(new_nb,rMap));
}  


# reorderVector
#
# input: data vector, offsets vector
# returns: vector of same length as input data vector
#          reordered according to offsets
#
reorderVector = function(x, rMap) {
  if (!is.vector(x)) return(NULL);
  N = length(x);              
  result = vector("numeric", length=N);
  for (i in 1:N) {
    result[i]= x[rMap[i]];
  }
  return(result);
}


# reorderMatrix
#
# input: data matrix, offsets vector
# returns: matrix of same shape as input data matrix,
#          rows reordered according to offsets
#
reorderMatrix = function(x, rMap) {
  if (!is.matrix(x)) return(NULL);
  N = nrow(x);
  result = matrix(nrow=N, ncol=ncol(x));
  for (i in 1:N) {
    result[i,]= x[rMap[i],];
  }
  return(result);
}


# scale_nb_components
#
# input: nb_object
# returns: vector of per-component scaling factor (for BYM2 model)
# scaling factor for singletons is 0
#
scale_nb_components = function(x) {
  N = length(x);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = indexByComponent(comp_ids);

  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  scales = vector("numeric", length=num_comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    scales[i] = 0.0;
    if (N_subregions > 1) {
      # get adj matrix for this component
      drops = comp_ids != i;
      nb_tmp = droplinks(x, drops);      
      nb_graph = nb2subgraph(nb_tmp, i, comp_ids, offsets);
      adj.matrix = sparseMatrix( i=nb_graph$node1, j=nb_graph$node2, x=1, dims=c(N_subregions,N_subregions), symmetric=TRUE);
      # compute ICAR precision matrix
      Q =  Diagonal(N_subregions, rowSums(adj.matrix)) - adj.matrix;
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N_subregions) * max(diag(Q)) * sqrt(.Machine$double.eps)
      # Compute the diagonal elements of the covariance matrix subject to the 
      # constraint that the entries of the ICAR sum to zero.
      Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,N_subregions),e=0))
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor = exp(mean(log(diag(Q_inv))))
      scales[i] = scaling_factor;
    }
  }
  return(scales);
}

library(INLA)
get_scaling_factor = function(nbs) {
  #Build the adjacency matrix
  adj.matrix = sparseMatrix(i=nbs$node1,j=nbs$node2,x=1,symmetric=TRUE)
  #The ICAR precision matrix (note! This is singular)
  Q=  Diagonal(nbs$N, rowSums(adj.matrix)) - adj.matrix

  #Add a small jitter to the diagonal for numerical stability (optional but recommended)
  Q_pert = Q + Diagonal(nbs$N) * max(diag(Q)) * sqrt(.Machine$double.eps)

  # Compute the diagonal elements of the covariance matrix subject to the
  # constraint that the entries of the ICAR sum to zero.
  #See the function help for further details.
  Q_inv = inla.qinv(Q_pert, constr=list(A = matrix(1,1,nbs$N),e=0))

  #Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    return((mean(log(diag(Q_inv)))))
}

