#!/usr/bin/env julia

using SIUnits
using SIUnits.ShortUnits

const zero_speed = Vector([eps()m/s, eps()m/s, eps()m/s])
const zero_locat = Vector([eps()m, eps()m, eps()m])

const G = 6.673e-11(N*(m^2)/(kg^2))

#typealias Vector{T} Array{T,1}
#typealias Matrix{T} Array{T,2}
immutable Body
    mass #::kg
    rad  #::m
    loc  ::Vector
    vel  ::Vector
end
typealias Universe Set{Body}

function prt(A)
    println("--- : $A")
    A
end

function length(v::Vector)
    sqrt(sum(v .^ 2))
end

function distance(v::Vector, w::Vector)
    length(v .- w)
end

function distance(b::Body, c::Body)
    distance(b.loc, c.loc)
end

function acc_by_grav(b::Body)
    b.mass*G/b.rad^2
end

function acc_by_grav(onto::Body, from::Body)
    d = distance(onto, from)
    from.mass*G/d^2
end

function unit(v::Vector)
    w = map(float, v)
    w/length(w)
end

function direction(b::Body, c::Body)
    d = c.loc .- b.loc
    unit(d)
end

function Update(b::Body, u::Universe, time)
    acting = setdiff(u, Set([b]))
    accel  = sum(map((c -> direction(b, c) .* acc_by_grav(b, c)), acting))

    pos = b.loc + (b.vel .* time) + (accel .* time ^ 2)

    velocity = b.vel .+ accel * time
    Body(b.mass, b.rad, pos, velocity)
end

function Update(u::Universe, time=86400s)
    Universe(map((b -> Update(b, u, time)), u))
end

earth = Body(5.972e24kg, 6371km,   Vector([1e10km,1e10km,1e10km]), zero_speed)
sun   = Body(1.989e30kg, 695700km, zero_locat, zero_speed)
earth

acc_by_grav(earth)
acc_by_grav(sun, earth)

I = Universe([earth, sun])

function steps(U::Universe, s::Integer)
    (s < 1) ? U : steps(Update(U), s-1)
end

steps(I, 100)
