
# just an example in case needed

struct mynewdistribution{T<:Real} <: Distribution{Univariate,Continuous} 
    a::T
    b::T
    # constructor
    function mynewdistribution{T}(a::T, b::T; check_args=true) where {T<:Real} 
        check_args && Distributions.@check_args(mynewdistribution, a>0 && b>0)
        return new{T}(a,b)
    end
end

# test
mynewdistribution{Float64}(2.0, 3.0)


# constructor for implicit calls
function mynewdistribution(a::Float64, b::Float64; check_args=true)
    return mynewdistribution{Float64}(a,b,check_args=check_args)
end

# constructor for different types
mynewdistribution(a::Real, b::Real) = mynewdistribution(promote(a,b,)...)
mynewdistribution(a::Integer, b::Integer) = mynewdistribution(float(a),float(b)) 

# test
mynewdistribution(2, 3)


import Base.rand, StatsBase.params
import Random, Distributions, Statistics, StatsBase

# helper func
params(d::mynewdistribution) = (d.a, d.b)

# rand
function Base.rand(rnd::AbstractRNG, d::mynewdistribution)
    (a, b) = params(d) 
    u = rand(rng)  # uniform
    return( quantile of mynewdistribution = Normal?? )
end

# sampler
Distributions.sampler(rng::AbstractRNG, d::mynewdistribution) = Base.rand(rng::AbstractRNG, d::mynewdistribution)


# logpdf
function Distributions.pdf(d::mynewdistribution{T}, x::Real) where {T<:Real}
    (a, b) = params(d)
    if x<=0 
        return zero(T)
    elseif x >= 1
        return zero(T)
    else
        return( PDF  ) 
    end
end

Distributions.logpdf(d::mynewdistribution, x::Real) = log(pdf(d,x))


function Statistics.quantile(d::mynewdistribution{T}, x::Real) where T<:Real
    (a, b) = params(d)
    if x<=0 
        return zero(T)
    elseif x >= 1
        return zero(T)
    else
        return( quantile function  ) 
    end
end

function Base.minimum(d::mynewdistribution) 
    return(0)
end

function Base.maximum(d::mynewdistribution) 
    return(1)
end

function Distributions.insupport(d::mynewdistribution)
    insupport(d::mynewdistribution, x::Real) = zero(x) <= x <= one(x)
end

Bijectors.bijector(d::mynewdistribution) = Logit(0., 1.)


if some_condition_to_reject
    Turing.@addlogprob! -Inf
    # Exit the model evaluation early
    return
end

Turing.@addlogprob! logpdf(x, μ)




