

function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end


function scaling_factor_bym2( adjacency_mat )
    # re-scaling variance using Reibler's solution and 
    # Buerkner's implementation: https://codesti.com/issue/paul-buerkner/brms/1241)  
    # Compute the diagonal elements of the covariance matrix subject to the 
    # constraint that the entries of the ICAR sum to zero.
    # See the inla.qinv function help for further details.
    # Q_inv = inla.qinv(Q, constr=list(A = matrix(1,1,nbs$N),e=0))  # sum to zero constraint
    # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
    # scaling_factor = exp(mean(log(diag(Q_inv))))
    asum = vec( sum(adjacency_mat, dims=2)) 
    N = size(g)[1]
    asum = float(asum) + N .* max.(asum) .* sqrt(1e-15)  # small perturbation
    Q = Diagonal(asum) 
    S = inv(Q)
    A = ones(size(Q)[1])   # constraint (sum to zero)
    V = S * A
    S = S - V * inv(A' * V) * V'
    scale_factor = exp(mean(log.(diag(S))))
    return scale_factor
end


function scaling_factor_bym2_groups(node1, node2, groups) 
    ## calculate the scale factor for each of k connected group of nodes, using the scale_c function from M. Morris
    gr = unique( groups )
    n_groups = length(gr)
    scale_factor = ones(n_groups)
    for j in 1:n_groups 
      k = findall( x -> x==j, groups)
      if length(k) > 1 
        e = Edge.(node1[k], node2[k])
        g = Graph(e)
        adjacency_mat = adjacency_matrix(g)
        scale_factor[j] = scaling_factor_bym2( adjacency_mat )
      end
    end
    return scale_factor
end
  


@model function turing_car(D, W, X, log_offset, y,  nbeta=size(X)[2] )
    # base model .. slow
    alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
    tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
    beta ~ filldist( Normal(0.0, 1.0), nbeta );
 
    prec = tau .* (D - alpha .* W)
    
    if !isposdef(prec)
        # check postive definiteness
        # phi ~ MvNormal( zeros(N) ); 
        Turing.@addlogprob! -Inf
        return nothing
    end
 
    phi ~ MvNormal(Symmetric(inv(prec)) );  # mean zero
    lambda = exp.( X * beta .+ phi .+ log_offset )
    @. y ~ Poisson( lambda );

end

   
@model function turing_car_prec(D, W, X, log_offset, y,  nbeta=size(X)[2] )
    alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
    tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
    beta ~ filldist( Normal(0.0, 1.0), nbeta);
 
    prec = tau .* (D - alpha .* W)
    
    if !isposdef(prec)
        # check postive definiteness
        # phi ~ MvNormal( zeros(N) ); 
        Turing.@addlogprob! -Inf
        return nothing
    end

    phi ~ MvNormalCanon( Symmetric(prec) );  # mean zero .. no inverse
    lambda = exp.( X * beta .+ phi .+ log_offset )
    @. y ~ Poisson( lambda );

end





@model function turing_icar_direct_test( node1, node2, ysd  )
    # equivalent to Morris' "simple_iar' .. testing pairwise difference formulation
    # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)

    phi ~ filldist( Normal(0.0, ysd), N)   # 10 is std from data: std(y)=7.9 stan goes from U(-Inf,Inf) .. not sure why 
    dphi = phi[node1] - phi[node2]
    lp_phi =  -0.5 * dot( dphi, dphi )
    Turing.@addlogprob! lp_phi
    
    # soft sum-to-zero constraint on phi)
    # equivalent to mean(phi) ~ normal(0,0.001)
    sum_phi = sum(phi)
    sum_phi ~ Normal(0, 0.001 * N);  
  
    # no data likelihood -- just prior sampling  -- 
end

  
@model function turing_icar_direct_bym( X, log_offset, y, node1, node2,  nbeta=size(X)[2], ysd=std(y) )
    # BYM
    # alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
     # tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
     beta ~ filldist( Normal(0.0, 5.0), nbeta);
     theta ~ filldist( Normal(0.0, 1.0), N) # unstructured (heterogeneous effect)
     # phi ~ filldist( Laplace(0.0, ysd), N) # spatial effects: stan goes from -Inf to Inf .. 
     phi ~ filldist( Normal(0.0, ysd), N) # spatial effects: stan goes from -Inf to Inf .. 
 
     # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
     # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
     dphi = phi[node1] - phi[node2]
     lp_phi =  -0.5 * dot( dphi, dphi )
     Turing.@addlogprob! lp_phi
     
     # soft sum-to-zero constraint on phi)
     # equivalent to mean(phi) ~ normal(0, 0.001)
     sum_phi = sum(phi)
     sum_phi ~ Normal(0, 0.001 * N);  

     tau_theta ~ Gamma(3.2761, 1.0/1.81);  # Carlin WinBUGS priors
     tau_phi ~ Gamma(1.0, 1.0);            # Carlin WinBUGS priors

     sigma_theta = inv(sqrt(tau_theta));  # convert precision to sigma
     sigma_phi = inv(sqrt(tau_phi));      # convert precision to sigma

     lambda = exp.( X * beta .+ phi .* sigma_phi .+ theta .* sigma_theta .+ log_offset )
  
     @. y ~ Poisson( lambda );
end
 

@model function turing_icar_direct_bym2(X, log_offset, y, node1, node2, scaling_factor,  nbeta=size(X)[2], ysd=std(y) )
    # BYM2
    # alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
     # tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
     beta ~ filldist( Normal(0.0, 5.0), nbeta);
     theta ~ filldist( Normal(0.0, 1.0), N)  # unstructured (heterogeneous effect)
#     phi ~ filldist( Laplace(0.0, ysd), N) # spatial effects: stan goes from -Inf to Inf .. 
     phi ~ filldist( Normal(0.0, ysd), N) # spatial effects: stan goes from -Inf to Inf .. 
        
     # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
     # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
     dphi = phi[node1] - phi[node2]
     lp_phi =  -0.5 * dot( dphi, dphi )
     Turing.@addlogprob! lp_phi
     
     # soft sum-to-zero constraint on phi)
     # equivalent to mean(phi) ~ normal(0, 0.001)
     sum_phi = sum(phi)
     sum_phi ~ Normal(0, 0.001 * N);  

     sigma ~ truncated( Normal(0, 1.0), 0, Inf) ; 
     rho ~ Beta(0.5, 0.5);

     # variance of each component should be approximately equal to 1
     convolved_re =  sigma .*  ( sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi );
   
     lambda = exp.( X * beta .+  convolved_re .+ log_offset )
  
     @. y ~ Poisson( lambda );

    # to compute from posteriors
    #  real logit_rho = log(rho / (1.0 - rho));
    #  vector[N] eta = log_E + beta0 + x * betas + convolved_re * sigma; // co-variates
    #  vector[N] lambda = exp(eta);
end
 

@model function turing_icar_direct_bym2_groups(X, log_offset, y, node1, node2, scaling_factor, groups, gi, nbeta=size(X)[2], ysd=std(y) )
    # BYM2
    # alpha ~ Uniform(0.0, 1.0); # alpha = 0.9 ; alpha==1 for BYM / iCAR
     # tau ~ Gamma(2.0, 1.0/2.0);  # tau=0.9
     beta ~ filldist( Normal(0.0, 5.0), nbeta);
     theta ~ filldist( Normal(0.0, 1.0), N)  # unstructured (heterogeneous effect)
     phi ~ filldist(Normal(0.0, ysd), N) # spatial effects: stan goes from -Inf to Inf .. 
        
     # pairwise difference formulation ::  prior on phi on the unit scale with sd = 1
     # see (https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
     dphi = phi[node1] - phi[node2]
     lp_phi =  -0.5 * dot( dphi, dphi )
     Turing.@addlogprob! lp_phi
     
     sigma ~ truncated( Normal(0, 1.0), 0, Inf) ; 
     rho ~ Beta(0.5, 0.5);

     convolved_re = zeros(N)

     for j in 1:length(gi)
         ic = gi[j] 
        
         # soft sum-to-zero constraint on phi)
         # equivalent to mean(phi) ~ normal(0, 0.001)
         sum_phi = sum(phi[ic])
         sum_phi ~ Normal(0, 0.001 * N);  

         if  length(ic) == 1 
             convolved_re[ ic ] = sigma .* theta[ ic ];
         else  
             convolved_re[ ic ] = sigma .* ( sqrt.(1 .- rho) .* theta[ ic ]  +  sqrt(rho ./ scaling_factor[j] )  .* phi[ ic ] ) ;
         end 
     end
  
     # convolved_re =  sqrt.(1 .- rho) .* theta .+ sqrt.(rho ./ scaling_factor) .* phi;
   
     lambda = exp.( X * beta .+  convolved_re .+ log_offset )
   
     @. y ~ Poisson( lambda );
  
    # to compute from posteriors
    #  real logit_rho = log(rho / (1.0 - rho));
    #  vector[N] eta = log_E + beta0 + x * betas + convolved_re * sigma; // co-variates
    #  vector[N] lambda = exp(eta);
end
  

function nodes( adj )
    N_edges = Integer( length(adj) / 2 );
    node1 =  fill(0, N_edges); 
    node2 =  fill(0, N_edges); 

    i_adjacency = 0;
    i_edge = 0;
    for i in 1:N
    for j in 1:num[i]
        i_adjacency = i_adjacency + 1;
        if i < adj[i_adjacency]
            i_edge = i_edge + 1;
            node1[i_edge] = i;
            node2[i_edge] = adj[i_adjacency];
        end
    end
    end
    return node1, node2
end
