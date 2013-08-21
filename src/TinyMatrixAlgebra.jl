#
# - Tiny Matrix Algebra -
#
# An implementation of matrix algebra for 3x1, 1x3, and 3x3 matricies.
#
# NOTE: 1X3 matrix is not yet implemented
#
# ==============================================================================
#
# Copyright 2013 Benjamin Silbaugh (ben.silbaugh@gmail.com)
#
# Permission is granted to redistribute under the terms of the MIT license.
#
# ==============================================================================

module TinyMatrixAlgebra

import Base.getindex
import Base.setindex!
import Base.cross
import Base.dot
import Base.sin
import Base.cos
import Base.tan
import Base.asin
import Base.acos
import Base.atan
import Base.+
import Base.-
import Base.*

export Matrix3X1, Matrix1X3, Matrix3X3,
       diag3, id3, 
       zero3X1, zero1X3, zero3X3,
       rand3X1, rand1X3, rand3X3,
       array,
       length, size, getindex, setindex!,
       +, -, *,
       inner, outer, dot, cross,
       sin, cos, asin, acos, tan, atan

# ------------------------------------------------------------------------------
# Matrix types
# ------------------------------------------------------------------------------

# 3X1 matrix (column vector)
immutable Matrix3X1{T <: Number}
    elem_11 :: T
    elem_21 :: T
    elem_31 :: T
end

# 1X3 matrix (row vector)
immutable Matrix1X3{T <: Number}
    elem_11 :: T
    elem_12 :: T
    elem_13 :: T
end

# 3X3 matrix (square matrix)
immutable Matrix3X3{T <: Number}
    elem_11 :: T
    elem_12 :: T
    elem_13 :: T
    elem_21 :: T
    elem_22 :: T
    elem_23 :: T
    elem_31 :: T
    elem_32 :: T
    elem_33 :: T
end

# ------------------------------------------------------------------------------
# Constructors/Initializers
# ------------------------------------------------------------------------------

# Diagonal matricies

function diag3(c::Number)
    o = zero(c)
    Matrix3X3(c, o, o,
              o, c, o,
              o, o, c)
end

id3(T::DataType) = diag3(one(T))

# Constant matricies

Matrix3X1(c::Number) = Matrix3X1(c, c, c)

Matrix1X3(c::Number) = Matrix1X3(c, c, c)

Matrix3X3(c::Number) = Matrix3X3(c, c, c, c, c, c, c, c, c)

zero3X1(T::DataType) = Matrix3X1(zero(T))

zero1X3(T::DataType) = Matrix1X3(zero(T))

zero3X3(T::DataType) = Matrix3X3(zero(T))

# Random matricies (limited to Float type returned by rand)

rand3X1() = Matrix3X1(rand(), rand(), rand())

rand1X3() = Matrix1X3(rand(), rand(), rand())

rand3X3() = Matrix3X3(rand(), rand(), rand(),
                      rand(), rand(), rand(),
                      rand(), rand(), rand())

# ------------------------------------------------------------------------------
# Transposition
# ------------------------------------------------------------------------------

# (could this be replaced with a reinterpret cast?
function transpose{T <: Number}(u::Matrix3X1{T}) 
    Matrix1X3(u.elem_11, u.elem_21, u.elem_31)
end

# (could this be replaced with a reinterpret cast?
function transpose{T <: Number}(u::Matrix1X3{T})
    Matrix3X1(u.elem_11, u.elem_12, u.elem_13)
end

function transpose{T <: Number}(a::Matrix3X3{T})
    Matrix3X3(a.elem_11, a.elem_21, a.elem_31,
              a.elem_12, a.elem_22, a.elem_32,
              a.elem_13, a.elem_23, a.elem_33)
end

# ------------------------------------------------------------------------------
# Mappings to/from Julia built-in types
# ------------------------------------------------------------------------------

# Copies the element of a (tiny) matrix into an equivolent Julia array

function array{T <: Number}(a::Matrix3X3{T})
    b = Array(T, 3, 3)
    for j = 1:3
        for i = 1:3
            b[i,j] = a[i,j]
        end
    end
    b
end

function array{T <: Number}(a::Matrix3X1{T})
    b = Array(T, 3)
    for i = 1:3
        b[i] = a[i]
    end
    b
end

function array{T <: Number}(a::Matrix1X3{T})
    b = Array(T, 1, 3)
    for i = 1:3
        b[i] = a[i]
    end
    b
end

# ------------------------------------------------------------------------------
# Higher-order (helper) functions
# ------------------------------------------------------------------------------

# Apply a operator/function element-wise
#
# Given a unitary operator op and matrix a, construct a new matrix of same
# dimension as a by setting each element of the new matrix equal to the
# result of applying op to the corresponding element of a.
#
# Example: Suppose a is the (3X3) Matrix. 
#
#     [a11, a12, a13]
#     [a21, a22, a23]
#     [a31, a32, a33]
#
# Then
#
#    b = apply(sin, a)
#
# results in b having the elements
#
#    [sin(a11), sin(a12), sin(a13)]
#    [sin(a21), sin(a22), sin(a23)]
#    [sin(a31), sin(a32), sin(a33)]
#
# Likewise, given a binary operator op and a pair of matricies a and b of equal
# dimension, construct a new matrix of same dimension by settings its elements
# equal to the result of applying the binary operator to the corresponding
# elements of a and b.
#
# Example: Suppose a and b are 3X3 matricies, and f is a function of two
# variables (of same type as the elements of a and b). Then
#
#     c = apply(f, a, b)
#
# results in c having the elements
#
#     [f(a11, b11), f(a12, b12), f(a13, b13)]
#     [f(a21, b21), f(a22, b22), f(a23, b23)]
#     [f(a31, b31), f(a32, b32), f(a33, b33)]
#
# NOTE:
# 1. In the case of a unitary operator, this behaves like the map function.

function apply{T <: Number}(op, a::Matrix3X1{T})
    Matrix3X1(op(a.elem_11),
              op(a.elem_21),
              op(a.elem_31))
end

function apply{T <: Number}(op, a::Matrix1X3{T})
    Matrix1X3(op(a.elem_11), op(a.elem_12), op(a.elem_13))
end

function apply{T <: Number}(op, a::Matrix3X3{T})
    Matrix3X3(op(a.elem_11),op(a.elem_12),op(a.elem_13),
              op(a.elem_21),op(a.elem_22),op(a.elem_23),
              op(a.elem_31),op(a.elem_32),op(a.elem_33))
end

function apply{T <: Number}(op, a::Matrix3X1{T}, b::Matrix3X1{T})
    Matrix3X1{T}( op(a.elem_11, b.elem_11), 
                  op(a.elem_21, b.elem_21),
                  op(a.elem_31, b.elem_31))
end

function apply{T <: Number}(op, a::Matrix1X3{T}, b::Matrix1X3{T})
    Matrix1X3( op(a.elem_11, b.elem_11),
               op(a.elem_12, b.elem_12),
               op(a.elem_13, b.elem_13))
end

function apply{T <: Number}(op, a::Matrix3X3{T}, b::Matrix3X3{T})
    Matrix3X3{T}( # first row
               op(a.elem_11, b.elem_11), 
               op(a.elem_12, b.elem_12),
               op(a.elem_13, b.elem_13),
               # second row
               op(a.elem_21, b.elem_21), 
               op(a.elem_22, b.elem_22),
               op(a.elem_23, b.elem_23),
               # third row
               op(a.elem_31, b.elem_31), 
               op(a.elem_32, b.elem_32),
               op(a.elem_33, b.elem_33))
end

# ------------------------------------------------------------------------------
# Container function overloads
# ------------------------------------------------------------------------------

# length overloads
# This should return the total number of elements

length{T <: Number}(u::Matrix3X1{T}) = 3

length{T <: Number}(u::Matrix1X3{T}) = 3

length{T <: Number}(u::Matrix3X3{T}) = 9

# size overloads
# This should return the matrix dimensions

size{T <: Number}(u::Matrix3X1{T}) = (3,1)

size{T <: Number}(u::Matrix1X3{T}) = (1,3)

size{T <: Number}(u::Matrix3X3{T}) = (3,3)

# Container indexing
# These may be slow -- don't use for performance critical code

macro index_error()
    :(throw(KeyError()))
end

function getindex{T <: Number}(u::Matrix3X1{T}, i::Int)
    1 == i ? u.elem_11 : 
    2 == i ? u.elem_21 : 
    3 == i ? u.elem_31 : @index_error
end

function getindex{T <: Number}(u::Matrix1X3{T}, i::Int)
    1 == i ? u.elem_11 : 
    2 == i ? u.elem_12 : 
    3 == i ? u.elem_13 : @index_error
end

function getindex{T <: Number}(u::Matrix3X3{T}, i::Int, j::Int)
    if 1 == j 
        1 == i ? u.elem_11 : 
        2 == i ? u.elem_21 : 
        3 == i ? u.elem_31 : @index_error
    elseif 2 == j
        1 == i ? u.elem_12 : 
        2 == i ? u.elem_22 : 
        3 == i ? u.elem_32 : @index_error
    elseif 3 == j
        1 == i ? u.elem_13 : 
        2 == i ? u.elem_23 : 
        3 == i ? u.elem_33 : @index_error
    else
        @index_error
    end
end

function setindex!{T <: Number}(u::Matrix3X1{T}, val::T, i::Int)
    1 == i ? u.elem_11 = val :
    2 == i ? u.elem_21 = val :
    3 == i ? u.elem_31 = val : @index_error
    u
end

function setindex!{T <: Number}(u::Matrix1X3{T}, val::T, i::Int)
    1 == i ? u.elem_11 = val :
    2 == i ? u.elem_12 = val :
    3 == i ? u.elem_13 = val : @index_error
    u
end

function setindex!{T <: Number}(u::Matrix3X3{T}, val::T, i::Int, j::Int)
    if 1 == j
        1 == i ? u.elem_11 = val :
        2 == i ? u.elem_21 = val :
        3 == i ? u.elem_31 = val : @index_error
    elseif 2 == j
        1 == i ? u.elem_12 = val :
        2 == i ? u.elem_22 = val :
        3 == i ? u.elem_32 = val : @index_error
    elseif 3 == j
        1 == i ? u.elem_13 = val :
        2 == i ? u.elem_23 = val :
        3 == i ? u.elem_33 = val : @index_error
    else
        @index_error
    end
    u
end

# ------------------------------------------------------------------------------
# Arithmetic
# ------------------------------------------------------------------------------

# Scalar multiplication and division
for mtype = (:Matrix3X1, :Matrix1X3, :Matrix3X3)
    eval(quote
        *{T <: Number}(u::($mtype){T}, c::T) = apply(x -> x*c, u)
        *{T <: Number}(c::T, u::($mtype){T}) = u*c
        /{T <: Number}(u::($mtype){T}, c::T) = apply(x -> x/c, u)
    end)
end

# Addition and subtraction
for mtype = (:Matrix3X1, :Matrix1X3, :Matrix3X3)
    for op = (:+, :-)
        eval(quote
            ($op){T <: Number}(u::($mtype){T}, v::($mtype){T}) = 
                apply($op, u, v)
            end)
    end
end

# ------------------------------------------------------------------------------
# Matrix products
# ------------------------------------------------------------------------------

function *{T <: Number}(a::Matrix3X3{T}, b::Matrix3X3{T})
    # first row
    elem_11 = a.elem_11*b.elem_11 + a.elem_12*b.elem_21 + a.elem_13*b.elem_31
    elem_12 = a.elem_11*b.elem_12 + a.elem_12*b.elem_22 + a.elem_13*b.elem_32
    elem_13 = a.elem_11*b.elem_13 + a.elem_12*b.elem_23 + a.elem_13*b.elem_33
    # second row
    elem_21 = a.elem_21*b.elem_11 + a.elem_22*b.elem_21 + a.elem_23*b.elem_31
    elem_22 = a.elem_21*b.elem_12 + a.elem_22*b.elem_22 + a.elem_23*b.elem_32
    elem_23 = a.elem_21*b.elem_13 + a.elem_22*b.elem_23 + a.elem_23*b.elem_33
    # third row
    elem_31 = a.elem_31*b.elem_11 + a.elem_32*b.elem_21 + a.elem_33*b.elem_31
    elem_32 = a.elem_31*b.elem_12 + a.elem_32*b.elem_22 + a.elem_33*b.elem_32
    elem_33 = a.elem_31*b.elem_13 + a.elem_32*b.elem_23 + a.elem_33*b.elem_33
    # Form final result as matrix
    Matrix3X3{T}(elem_11, elem_12, elem_13, 
               elem_21, elem_22, elem_23, 
               elem_31, elem_32, elem_33)
end

function *{T <: Number}(a::Matrix3X3{T}, u::Matrix3X1{T})
    elem_11 = a.elem_11*u.elem_11 + a.elem_12*u.elem_21 + a.elem_13*u.elem_31
    elem_21 = a.elem_21*u.elem_11 + a.elem_22*u.elem_21 + a.elem_23*u.elem_31
    elem_31 = a.elem_31*u.elem_11 + a.elem_32*u.elem_21 + a.elem_33*u.elem_31
    Matrix3X1{T}(elem_11, elem_21, elem_31)
end

function *{T <: Number}(u::Matrix1X3{T}, a::Matrix3X3{T})
    elem_11 = u.elem_11*a.elem_11 + u.elem_12*a.elem_21 + u.elem_13*a.elem_31
    elem_12 = u.elem_11*a.elem_12 + u.elem_12*a.elem_22 + u.elem_13*a.elem_32
    elem_13 = u.elem_11*a.elem_13 + u.elem_12*a.elem_23 + u.elem_13*a.elem_33
    Matrix1X3(elem_11, elem_12, elem_13)
end

function *{T <: Number}(u::Matrix3X1{T}, v::Matrix1X3{T})
    # first row
    elem_11 = u.elem_11*v.elem_11
    elem_12 = u.elem_11*v.elem_12
    elem_13 = u.elem_11*v.elem_13
    # second row
    elem_21 = u.elem_21*v.elem_11
    elem_22 = u.elem_21*v.elem_12
    elem_23 = u.elem_21*v.elem_13
    # third row
    elem_31 = u.elem_31*v.elem_11
    elem_32 = u.elem_31*v.elem_12
    elem_33 = u.elem_31*v.elem_13
    # result
    Matrix3X3(elem_11, elem_12, elem_13,
              elem_21, elem_22, elem_23,
              elem_31, elem_32, elem_33)
end

function *{T <: Number}(u::Matrix1X3{T}, v::Matrix3X1{T})
    u.elem_11*v.elem_11 + u.elem_12*v.elem_21 + u.elem_13*v.elem_31
end

# Inner (dot) product

function inner{T <: Number}(u::Matrix3X1{T}, v::Matrix3X1{T})
    u.elem_11*v.elem_11 + u.elem_21*v.elem_21 + u.elem_31*v.elem_31
end

dot{T <: Number}(u::Matrix3X1{T}, v::Matrix3X1{T}) = inner(u,v)

# Outer product

function outer{T <: Number}(u::Matrix3X1{T}, v::Matrix3X1{T})
    # first row
    elem_11 = u.elem_11*v.elem_11
    elem_12 = u.elem_11*v.elem_21
    elem_13 = u.elem_11*v.elem_31
    # second row
    elem_21 = u.elem_21*v.elem_11
    elem_22 = u.elem_21*v.elem_21
    elem_23 = u.elem_21*v.elem_31
    # third row
    elem_31 = u.elem_31*v.elem_11
    elem_32 = u.elem_31*v.elem_21
    elem_33 = u.elem_31*v.elem_31
    # Form resulting matrix
    Matrix3X3{T}(elem_11, elem_12, elem_13,
               elem_21, elem_22, elem_23,
               elem_31, elem_32, elem_33)
end

# Cross product
# (We're only dealing with 3D, so this makes sense)

function cross{T <: Number}(a::Matrix3X1{T}, b::Matrix3X1{T})
    elem_11 =  a.elem_21*b.elem_31 - a.elem_31*b.elem_21
    elem_21 = -a.elem_11*b.elem_31 + a.elem_31*b.elem_11
    elem_31 =  a.elem_11*b.elem_21 - a.elem_21*b.elem_11
    Matrix3X1{T}(elem_11, elem_21, elem_31)
end

# ------------------------------------------------------------------------------
# Trigonometric functions
# ------------------------------------------------------------------------------

for mtype = (:Matrix3X1, :Matrix1X3, :Matrix3X3)
    for fun = (:sin, :cos, :tan, :asin, :acos, :atan)
        eval(quote
            ($fun){T}(a::($mtype){T}) = apply(($fun), a)
        end)
    end
end

end #module
