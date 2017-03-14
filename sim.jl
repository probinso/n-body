#!/usr/bin/env julia

using SIUnits
using SIUnits.ShortUnits

using Plots

const zero_vel = Vector([eps()m/s, eps()m/s, eps()m/s])
const zero_pos = Vector([eps()m, eps()m, eps()m])
const G = 6.673e-11(N*(m^2)/(kg^2))

immutable Body
    mass #::kg
    rad  #::m
    pos  ::Vector
    vel  ::Vector
end

const zero_body = Body(eps()kg, eps()km, zero_pos, zero_vel)

typealias Universe Set{Body}

function prt(A)
    println("--- : $A")
    A
end

function vlength(v::Vector)
    sqrt(sum(v .^ 2))
end

function distance(v::Vector, w::Vector)
    vlength(v .- w)
end

function distance(b::Body, c::Body)
    distance(b.pos, c.pos)
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
    w/vlength(w)
end

function direction(b::Body, c::Body)
    d = c.pos .- b.pos
    unit(d)
end

function Update(b::Body, u::Universe, time)
    acting = setdiff(u, Set([b]))

    acc = sum(map((c -> direction(b, c) .* acc_by_grav(b, c)), acting))
    vel = b.vel .+ acc * time
    pos = b.pos .+ (b.vel .* time) .+ (acc .* time ^ 2)

    Body(b.mass, b.rad, pos, vel)
end

function Join(a::Body, b::Body)
    Body(
         a.mass + b.mass,
         sum([a.rad, b.rad] .^ 3)^(1/3),
         (a.pos + b.pos)/2,
         (a.mass*a.vel + b.mass*b.vel)/(a.mass + b.mass)
         )
end

function collisions(u::Universe)
    group_collide = (b, u) -> filter(c -> distance(b, c) < max(b.rad, c.rad), u)

    v = Universe()
    while !isempty(u)
        b = pop!(u)
        group = group_collide(b, u)
        for c in group
            delete!(u, c)
        end
        group = union(group, [b])
        v = union(v, [reduce(Join, group)])
    end
    Universe(v)
end

function Update(u::Universe, time=3600s)
    v = collisions(u)
    Universe(map(b -> Update(b, v, time), v))
end

sun =
    Body(
         1.989e30kg,
         695700km,
         zero_pos,
         zero_vel
         )

mercury =
    Body(
         3.285e23kg,
         2440km,
         zero_pos,
         zero_vel
         )

venus =
    Body(
         4.867e24kg,
         6052km,
         zero_pos,
         zero_vel
         )

earth =
    Body(
         5.972e24kg,
         6371km,
         Vector([1e10km,1e10km,1e10km]),
         zero_vel
         )

mars =
    Body(
         6.39e23kg,
         3390km,
         zero_pos,
         zero_vel
         )

jupiter =
    Body(
         1.898e27kg,
         69911km,
         zero_pos,
         zero_vel
         )

saturn  =
    Body(
         5.683e26kg,
         58232km,
         zero_pos,
         zero_vel
         )

uranus  = zero_body
neptune = zero_body



acc_by_grav(earth)
acc_by_grav(sun, earth)

I = Universe([sun, mercury, venus, earth, mars, jupiter, saturn])

function steps(U::Universe, s::Integer)
    (s < 1) ? U : steps(Update(U), s-1)
end

function get_positions(u::Universe)
    lround = f -> round(f, 3)
    map(b -> map(lround, map(float, b.pos)), u)
end



steps(I, 10)

# using Distributions
