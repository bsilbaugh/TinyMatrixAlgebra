#
#  - Tiny Matrix Algebra -
#
# Correctness tests
#
# Copyright 2013 Benjamin Silbaugh (ben.silbaugh@gmail.com)
#
# Permission is granted to redistributed under the terms of the MIT license
#
# ==============================================================================

include("../src/TinyMatrixAlgebra.jl")
include("testing.jl")

using TinyMatrixAlgebra

typealias Float Float64
typealias Matrix Matrix3X3{Float}
typealias Vector Matrix3X1{Float}
typealias Covector Matrix1X3{Float}

function length_vector()
    u = rand3X1()
    @check_eq 3 length(u)
end

function length_covector()
    u = rand1X3()
    @check_eq 3 length(u)
end

function length_matrix()
    u = rand3X3()
    @check_eq 9 length(u)
end

function size_vector()
    u = rand3X1()
    @check_eq (3,1) size(u)
end

function size_covector()
    u = rand1X3()
    @check_eq (1,3) size(u)
end

function size_matrix()
    u = rand3X3()
    @check_eq (3,3) size(u)
end

function getindex_vector()
    const elem_11 = rand()
    const elem_21 = rand()
    const elem_31 = rand()
    u = Vector(elem_11, elem_21, elem_31)
    @check_eq elem_11 u[1]
    @check_eq elem_21 u[2]
    @check_eq elem_31 u[3]
end

function getindex_covector()
    const elem_11 = rand()
    const elem_12 = rand()
    const elem_13 = rand()
    u = Covector(elem_11, elem_12, elem_13)
    @check_eq elem_11 u[1]
    @check_eq elem_12 u[2]
    @check_eq elem_13 u[3]
end

function getindex_matrix()
    const elem_11 = rand()
    const elem_12 = rand()
    const elem_13 = rand()
    const elem_21 = rand()
    const elem_22 = rand()
    const elem_23 = rand()
    const elem_31 = rand()
    const elem_32 = rand()
    const elem_33 = rand()
    a = Matrix(elem_11, elem_12, elem_13,
               elem_21, elem_22, elem_23,
               elem_31, elem_32, elem_33)

    @check_eq elem_11 a[1,1]
    @check_eq elem_12 a[1,2]
    @check_eq elem_13 a[1,3]

    @check_eq elem_21 a[2,1]
    @check_eq elem_22 a[2,2]
    @check_eq elem_23 a[2,3]

    @check_eq elem_31  a[3,1]
    @check_eq elem_32  a[3,2]
    @check_eq elem_33  a[3,3]

end

function setindex_vector()
    u = zero3X1(Float)
    u[1] = 0.1
    u[2] = 0.2
    u[3] = 0.3
    @check_eq 0.1  u[1]
    @check_eq 0.2  u[2]
    @check_eq 0.3  u[3]
end

function setindex_covector()
    u = zero1X3(Float)
    u[1] = 0.1
    u[2] = 0.2
    u[3] = 0.3
    @check_eq 0.1  u[1]
    @check_eq 0.2  u[2]
    @check_eq 0.3  u[3]
end

function setindex_matrix()
    a = zero3X3(Float)
    a[1,1] = 0.11
    a[1,2] = 0.12
    a[1,3] = 0.13
    a[2,1] = 0.21
    a[2,2] = 0.22
    a[2,3] = 0.23
    a[3,1] = 0.31
    a[3,2] = 0.32
    a[3,3] = 0.33
    @check_eq a[1,1]  0.11
    @check_eq a[1,2]  0.12
    @check_eq a[1,3]  0.13
    @check_eq a[2,1]  0.21
    @check_eq a[2,2]  0.22
    @check_eq a[2,3]  0.23
    @check_eq a[3,1]  0.31
    @check_eq a[3,2]  0.32
    @check_eq a[3,3]  0.33
end

function scalar_multiply_vector()
    a = rand3X1()
    a0 = array(a)
    c = rand()
    b = c*a
    b0 = c*a0
    @check_eq_vector b0 b
end

function scalar_multiply_covector()
    a = rand1X3()
    a0 = array(a)
    c = rand()
    b = c*a
    b0 = c*a0
    @check_eq_vector b0 b
end

function scalar_multiply_matrix()
    a = rand3X3()
    a0 = array(a)
    c = rand()
    b = c*a
    b0 = c*a0
    @check_eq_matrix b0 b
end

function scalar_division_vector()
    u = rand3X1()
    u0 = array(u)
    c = 0.31432 # a "random" non-zero number
    v = u/c
    v0 = u0/c
    @check_eq_vector v0 v
end

function scalar_division_covector()
    u = rand3X1()
    u0 = array(u)
    c = 0.31432 # a "random" non-zero number
    v = u/c
    v0 = u0/c
    @check_eq_vector v0 v
end

function scalar_division_matrix()
    a = rand3X3()
    a0 = array(a)
    c = 0.43135
    b = a/c
    b0 = a0/c
    @check_eq_matrix b0 b
end

function addition_vector()
    a = rand3X1()
    b = rand3X1()
    a0 = array(a)
    b0 = array(b)
    c = a + b
    c0 = a0 + b0
    @check_eq_vector c0 c
end

function addition_covector()
    a = rand1X3()
    b = rand1X3()
    a0 = array(a)
    b0 = array(b)
    c = a + b
    c0 = a0 + b0
    @check_eq_vector c0 c
end

function subtraction_vector()
    a = rand3X1()
    b = rand3X1()
    a0 = array(a)
    b0 = array(b)
    c = a - b
    c0 = a0 - b0
    @check_eq_vector c0 c
end

function subtraction_covector()
    a = rand1X3()
    b = rand1X3()
    a0 = array(a)
    b0 = array(b)
    c = a - b
    c0 = a0 - b0
    @check_eq_vector c0 c
end

function addition_matrix()
    a = rand3X3()
    b = rand3X3()
    a0 = array(a)
    b0 = array(b)
    c = a + b
    c0 = a0 + b0
    @check_eq_matrix c0 c
end

function subtraction_matrix()
    a = rand3X3()
    b = rand3X3()
    a0 = array(a)
    b0 = array(b)
    c = a - b
    c0 = a0 - b0
    @check_eq_matrix c0 c
end

function multiply_matrix_matrix()
    a = rand3X3()
    b = rand3X3()
    a0 = array(a)
    b0 = array(b)
    c = a*b
    c0 = a0*b0
    @check_eq_matrix c0 c
end

function multiply_matrix_vector()
    a = rand3X3()
    v = rand3X1()
    a0 = array(a)
    v0 = array(v)
    c = a*v
    c0 = a0*v0
    @check_eq_vector c0 c
end

function multiply_covector_matrix()
    a = rand3X3()
    v = rand1X3()
    a0 = array(a)
    v0 = array(v)
    c = v*a
    c0 = v0*a0
    @check_eq_vector c0 c
end

function multiply_covector_vector()
    v = rand1X3()
    u = rand3X1()
    v0 = array(v)
    u0 = array(u)
    c = v*u
    c0 = v0*u0
    @check_eq c0 c
end

function multiply_vector_covector()
    v = rand3X1()
    u = rand1X3()
    v0 = array(v)
    u0 = array(u)
    a = v*u
    a0 = v0*u0
    @check_eq_matrix a0 a
end

function inner_vector_vector()
    u = rand3X1()
    v = rand3X1()
    u0 = array(u)
    v0 = array(v)
    c = inner(u,v)
    c0 = dot(u0,v0)
    @check_eq c0 c
end

function outer_vector_vector()
    u = rand3X1()
    v = rand3X1()
    u0 = Array(Float64, 3, 1)
    v0 = Array(Float64, 1, 3)
    u0[:,1] = {u[1], u[2], u[3]}
    v0[1,:] = {v[1], v[2], v[3]}
    c = outer(u, v)
    c0 = u0*v0
    @check_eq_matrix c0 c
end

function cross_vector_vector()
    u = rand3X1()
    v = rand3X1()
    w = cross(u, v)
    u0 = array(u)
    v0 = array(v)
    w0 = cross(u0, v0)
    @check_eq_vector w0 w
end

@run_test length_vector
@run_test length_covector
@run_test length_matrix
@run_test size_vector
@run_test size_covector
@run_test size_matrix
@run_test getindex_vector
@run_test getindex_covector
@run_test getindex_matrix
@run_test setindex_vector
@run_test setindex_covector
@run_test setindex_matrix
@run_test scalar_multiply_vector
@run_test scalar_multiply_covector
@run_test scalar_multiply_matrix
@run_test scalar_division_vector
@run_test scalar_division_covector
@run_test scalar_division_matrix
@run_test addition_vector
@run_test addition_covector
@run_test subtraction_vector
@run_test subtraction_covector
@run_test addition_matrix
@run_test subtraction_matrix
@run_test multiply_matrix_matrix
@run_test multiply_matrix_vector
@run_test multiply_covector_matrix
@run_test multiply_covector_vector
@run_test multiply_vector_covector
@run_test inner_vector_vector
@run_test outer_vector_vector
@run_test cross_vector_vector
