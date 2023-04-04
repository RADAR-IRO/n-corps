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



function test(n::Int64)
    x = [Body(1, 0, 0, 1,1)]
    for t ∈ 2:n
        push!(x, predictorStep(x[t-1], 0.01))
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



