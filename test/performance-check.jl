#
# - Tiny Matrix Algebra -
#
# Performance tests
#
# Copyright 2013 Benjamin Silbaugh (ben.silbaugh@gmail.com)
#
# Permission is granted to redistribute under the terms of the MIT license
#
# ==============================================================================

require("../src/TinyMatrixAlgebra.jl")

using TinyMatrixAlgebra

import Base.cross
import TinyMatrixAlgebra.cross
import Base.dot
import TinyMatrixAlgebra.dot

function init3X3(n::Int)
    for i = 1:n
        rand3X3()
    end
end

function initArray(n::Int)
    for i = 1:n
        rand(3,3)
    end
end

function matmul(a, u, n::Int)
    for i = 1:n
        v = a*u
    end
end

function cross(u, v, n::Int)
    for i = 1:n
        w = cross(u, v)
    end
end

function dot(u, v, n::Int)
    for i = 1:n
        c = dot(u, v)
    end
end

n = 10000

a = rand3X3()
b = rand3X3()
u = rand3X1()
v = rand3X1()

a0 = array(a)
b0 = array(b)
u0 = array(u)
v0 = array(v)

println("3X3 initialization")
init3X3(n)
initArray(n)
@time init3X3(n)
@time initArray(n)

println("matrix-vector multiplication:")
matmul(a, u, n)
matmul(a0, u0, n)
@time matmul(a, u, n)
@time matmul(a0, u0, n)

println("matrix-matrix multiplication:")
matmul(a, b, n)
matmul(a0, b0, n)
@time matmul(a, b, n)
@time matmul(a0, b0, n)

println("cross product:")
cross(u, v, n)
cross(u0, v0, n)
@time cross(u, v, n)
@time cross(u0, v0, n)

println("dot product:")
dot(u, v, n)
dot(u0, v0, n)
@time dot(u, v, n)
@time dot(u0, v0, n)
