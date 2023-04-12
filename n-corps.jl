#!/bin/env julia

using StaticArrays #https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/
using Plots

# Single step of integrator: predictor then corrector

# Predictor
# first part: calculate x̃

# x has 4 parts: (r, θ, p, l)
struct Body <: FieldVector{5, Float64}
    r::Float64
    θ::Float64
    p::Float64
    l::Float64
    m::Float64
end

struct Ξ <: FieldVector{3, Float64}
    ξ1::Float64
    ξ2::Float64
    ξ3::Float64
end

# f also has 4 parts
function f(x::Body)
    x̃r = x.p/x.m #x.vr = x.p/x.m
    x̃θ = x.l/(x.r^2*x.m)
    x̃p = x.l^2/(x.r^3*x.m) - 1.5/x.r^2
    return Body(x̃r, x̃θ, x̃p, 0, 0)
end

# x̃ = x_o + τ * f(x_0, t)
function predictorStep(x_0::Body, τ::Float64)
    return x_0 + τ * f(x_0)
end

function T(x::Body)
    ξ1 = -1.5/x.r
    ξ2 = x.p^2/(2*x.m) + x.l^2/(2*x.m*x.r^2)
    ξ3 = x.l
    return Ξ(ξ1, ξ2, ξ3)
end

function Tprime(x::Body)
    dξ_1dr = 1.5/x.r^2
    dξ_2dr = -x.l^2/x.m * 1/x.r^3
    dξ_2dp = x.p/x.m
    dξ_2dl = x.l/(x.m * x.r^2)
    dξ_3dl = 1
    return [ dξ_1dr 0 0 0 0 ; dξ_2dr 0 dξ_2dp dξ_2dl 0 ; 0 0 0 dξ_3dl 0 ]
end

function solveForTheta(θ_0::Float64, l, m, r_0, r, k, v)
    A = l^2/(m*r_0) - k
    fθ_0 = A*v*cos(θ_0) + A*l/(m*r)*sin(θ_0) + k*v
    fprimeθ_0 = (A+k)*sin(θ_0) - k*v*sin(θ_0) + l^3/(m^2*r*r_0)*cos(θ_0) - k*l/(m*r)*cos(θ_0)
    θ_1 = θ_0 - fθ_0/fprimeθ_0
    if θ_1 - θ_0 <= 0
        return θ_1
    else 
        return solveForTheta(θ_1, l, m, r_0, r, k, v)
    end
end

function Tinverse(ξ::Ξ,x̃p,m,x̃θ,r_0,k)
    r = -k/ξ.ξ1
    l = ξ.ξ3
    p = sign(x̃p) * sqrt(2*m*ξ.ξ2-l^2/r^2)
    θ = solveForTheta(x̃θ, l, m, r_0, r, k, p/m)
    return Body(r, θ, p, l, m)
end

function integrator(x_0::Body, τ::Float64)
    x̃ = predictorStep(x_0, τ)
    ξ_0 = T(x_0)
    ξ_1 = ξ_0 + 0.5 * τ * (Tprime(x_0)*f(x_0)+Tprime(x̃)*f(x̃))
    return Tinverse(ξ_1, x̃.p, x̃.m, x̃.θ, x_0.r, 1.5)
end


function test(n::Int64)
    x = [Body(1, 0, 0, 1,1)]
    for t ∈ 2:n
        push!(x, integrator(x[t-1], 0.01))
    end
    return x
end

test(2)
bs = test(1000)

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

t = 1:1000
x = (v -> bs[v].r*cos(bs[v].θ)).(t)
y = (v -> bs[v].r*sin(bs[v].θ)).(t)

anim = @animate for i ∈ 1:1000
    circleplot(x, y, i)
end
gif(anim, "anim_fps60.gif", fps = 60)



