
# S ~ MVN(μ, vΣ)
# vΣ = v(I - γ C)^(-1) M
# I = N x N identity matrix
# M = N x N diagonal matrix, with elements Mii proportional to the conditional variance of Si | Sj
# C = N x N weight matrix, with elements Cij reflecting spatial association between areas i and j
# γ = controls overall strength of spatial dependence

# Si | S-i ~ Normal (μ i + Σj γ Cij (sj -μi), φ Mii ), 
# (Σj denotes summation)
# (-i all but i)

# in bym (improper CAR :: (the overall mean of the Si is not defined) ):
# Si | S-i ~ Normal ( S.bari , v / ni )
# S.bari = Σj in δi Sj / ni
# δi denotes the set of labels of the "neighbours" of area i.

# Proper CAR:
# The choice Cij = 1/ni if areas i and j are adjacent and Cij = 0 otherwise (with Cii also set to 0) , and 
# Mii = 1/ni, where ni is the number of areas which are adjacent to area i is still valid for the proper CAR

# In the context of disease mapping Cressie and Chan (1989) and Stern and Cressie (1999) choose
#  an alternative parameterisation:

#    Mii = 1/Ei (the inverse of the expected count or population size in area i)
#    Cij = (Ej / Ei)^(1/2) for neighbouring areas i, j and 0 otherwise

# MCAR:
# a multivariate p-dimensional vector of spatially correlated Gaussian data or
#  random effects in each area, Si = (S1i , S2i,....., Spi )', i=1,..., N.
#  Si | (S1(-i), S2(-i)) ~ Bivariate Normal (S.bari , V / ni ),
#  (here (S1(-i), S2(-i)) denotes the elements of the 2 x N matrix S excluding the ith area (column))
#  where S.bari = (S.bari1 ,S.bari2) with S.barip = Σ j in δi Sjp / ni and, as in 
#  the univariate case, δi and ni denote the set of labels of the "neighbours" of area i and 
#  the number of neighbours, respectively. V is a 2 x 2 covariance matix with diagonal 
#  elements v11 and v22 representing the conditional variances of S1 and S2 respectively, 
#  and off-diagonal element v12 representing the (conditional) within-area covariance between S1 and S2. 


# ϕi∣ϕj,j≠i∼N(α∑j=1nbijϕj,τ−1i)
# where τi is a spatially varying precision parameter, and bii=0.

# @model function bym()
#     phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
#     beta ~ normal(0, 1);
#     tau ~ gamma(2, 2);
#     y ~ poisson_log(X * beta + phi + log_offset);
# end
 
 
using Turing 
using PDMats
using LinearAlgebra 
using StatsModels
using DataFrames
using SparseArrays
using Graph

# stan implementation:
# https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html#auto-regressive-models-for-areal-data

# data source:  
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html

# y: the observed lip cancer case counts on a per-county basis
# x: an area-specific continuous covariate that represents the proportion of the population employed in agriculture, fishing, or forestry (AFF)
# E: the expected number of cases, used as an offset,
# adj: a list of region ids for adjacent regions
# num: a list of the number of neighbors for each region

N = 56

y   = [ 9, 39, 11,  9, 15,  8, 26,  7,  6, 20,
13,  5,  3,  8, 17,  9,  2,  7,  9,  7,
16, 31, 11,  7, 19, 15,  7, 10, 16, 11,
5,  3,  7,  8, 11,  9, 11,  8,  6,  4,
10,  8,  2,  6, 19,  3,  2,  3, 28,  6,
1,  1,  1,  1,  0,  0]

E = [1.4, 8.7, 3.0, 2.5, 4.3, 2.4, 8.1, 2.3, 2.0, 6.6,
4.4, 1.8, 1.1, 3.3, 7.8, 4.6, 1.1, 4.2, 5.5, 4.4,
10.5,22.7, 8.8, 5.6,15.5,12.5, 6.0, 9.0,14.4,10.2,
4.8, 2.9, 7.0, 8.5,12.3,10.1,12.7, 9.4, 7.2, 5.3,
18.8,15.8, 4.3,14.6,50.7, 8.2, 5.6, 9.3,88.7,19.6,
3.4, 3.6, 5.7, 7.0, 4.2, 1.8]

x = [16,16,10,24,10,24,10, 7, 7,16, 
7,16,10,24, 7,16,10, 7, 7,10,
7,16,10, 7, 1, 1, 7, 7,10,10,
7,24,10, 7, 7, 0,10, 1,16, 0,
1,16,16, 0, 1, 7, 1, 1, 0, 1,
1, 0, 1, 1,16,10]

adj = [ 5, 9,11,19,
7,10,
6,12,
18,20,28,
1,11,12,13,19,
3, 8,
2,10,13,16,17,
6,
1,11,17,19,23,29,
2, 7,16,22,
1, 5, 9,12,
3, 5,11,
5, 7,17,19,
31,32,35,
25,29,50,
7,10,17,21,22,29,
7, 9,13,16,19,29,
4,20,28,33,55,56,
1, 5, 9,13,17,
4,18,55,
16,29,50,
10,16,
9,29,34,36,37,39,
27,30,31,44,47,48,55,56,
15,26,29,
25,29,42,43,
24,31,32,55,
4,18,33,45,
9,15,16,17,21,23,25,26,34,43,50,
24,38,42,44,45,56,
14,24,27,32,35,46,47,
14,27,31,35,
18,28,45,56,
23,29,39,40,42,43,51,52,54,
14,31,32,37,46,
23,37,39,41,
23,35,36,41,46,
30,42,44,49,51,54,
23,34,36,40,41,
34,39,41,49,52,
36,37,39,40,46,49,53,
26,30,34,38,43,51,
26,29,34,42,
24,30,38,48,49,
28,30,33,56,
31,35,37,41,47,53,
24,31,46,48,49,53,
24,44,47,49,
38,40,41,44,47,48,52,53,54,
15,21,29,
34,38,42,54,
34,40,49,54,
41,46,47,49,
34,38,49,51,52,
18,20,24,27,56,
18,24,30,33,45,55]

weights = [ 1, 1,1,1,
1,1,
1,1,
1,1,1,
1,1,1,1,1,
1, 1,
1,1,1,1,1,
1,
1,1,1,1,1,1,
1, 1,1,1,
1, 1, 1,1,
1, 1,1,
1, 1,1,1,
1,1,1,
1,1,1,
1,1,1,1,1,1,
1, 1,1,1,1,1,
1,1,1,1,1,1,
1, 1, 1,1,1,
1,1,1,
1,1,1,
1,1,
1,1,1,1,1,1,
1,1,1,1,1,1,1,1,
1,1,1,
1,1,1,1,
1,1,1,1,
1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,
1,1,1,1,1,1,1,
1,1,1,1,
1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,1,1,
1,1,1,1,
1,1,1,1,1,
1,1,1,1,1,1,
1,1,1,1,1,
1,1,1,1,1,
1,1,1,1,1,1,1,
1,1,1,1,1,1,
1,1,1,1,
1,1,1,1,1,
1,1,1,1,
1,1,1,1,1,1,
1,1,1,1,1,1,
1,1,1,1,
1,1,1,1,1,1,1,1,1,
1,1,1,
1,1,1,1,
1,1,1,1,
1,1,1,1,
1,1,1,1,1,
1,1,1,1,1,
1,1,1,1,1,1]

             
num = [4, 2, 2, 3, 5, 2, 5, 1,  6, 
4, 4, 3, 4, 3, 3, 6, 6, 6 ,5, 
3, 3, 2, 6, 8, 3, 4, 4, 4,11,  
6, 7, 4, 4, 9, 5, 4, 5, 6, 5, 
5, 7, 6, 4, 5, 4, 6, 6, 4, 9, 
3, 4, 4, 4, 5, 5, 6]

N_edges = Integer( length(adj) / 2 );
node1 =  fill(0, N_edges); 
node2 =  fill(0, N_edges); 

iAdj = 0;
iEdge = 0;
for i in 1:N
for j in 1:num[i]
    iAdj = iAdj + 1;
    if i < adj[iAdj]
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adj[iAdj];
    end
end
end

e = Edge.(node1, node2)
g = Graph(e)
W = adjacency_matrix(g)

D = diagm(vec( sum(W, dims=2) ))

# x = collect(x)
x_scaled = (x .- mean(x)) ./ std(x)
X = DataFrame( Intercept=ones( N ), x=x_scaled )
X = Matrix(X)

log_offset = log.(E)

