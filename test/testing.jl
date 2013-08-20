#
# Testing
#
# A collection of ad-hoc utilties for testing julia code.
#
# Copyright 2013 Benjamin Silbaugh (ben.silbaugh@gmail.com)
#
# Permission is granted to redistribute under the terms of the MIT license.
#
# ==============================================================================

# NOTE: Some of these may yield false positives when applied to floating
# point numbers. Ideally, need to implement a @check_eq_close macro.

macro check_eq(a, b)
    checking_msg = string("check ", a, " == ", b)
    quote
    $a == $b ? println($checking_msg, " passed") : 
               println($checking_msg, " <----- FAILED! ", $a, "!=", $b)
    end
end

macro check_eq_matrix(a, b)
    quote
    @check_eq ($a)[1,1]  ($b)[1,1]
    @check_eq ($a)[1,2]  ($b)[1,2]
    @check_eq ($a)[1,3]  ($b)[1,3]
    @check_eq ($a)[2,1]  ($b)[2,1]
    @check_eq ($a)[2,2]  ($b)[2,2]
    @check_eq ($a)[2,3]  ($b)[2,3]
    @check_eq ($a)[3,1]  ($b)[3,1]
    @check_eq ($a)[3,2]  ($b)[3,2]
    @check_eq ($a)[3,3]  ($b)[3,3]
    end
end

macro check_eq_vector(u, v)
    quote
        @check_eq ($u)[1]  ($v)[1]
        @check_eq ($u)[2]  ($v)[2]
        @check_eq ($u)[3]  ($v)[3]
    end
end

macro run_test(test_name)
    quote
        println("=== running ", $test_name, " ===")
        $test_name()
    end
end
