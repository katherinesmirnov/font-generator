"""
    Creates a B-spline of knots "t" (list), centered at knot "i", and of degree "k"-1
    Note that indices start at 1, so i >= 1.
"""
function bspline(t, i, k)
    if k == 0
        f(x) = t[i] <= x < t[i+1] ? 1 : 0
        return f
    else
        f1(x) = (x-t[i])/(t[i+k]-t[i]) * bspline(t, i, k-1)(x)
        f2(x) = (t[i+k+1]-x)/(t[i+k+1]-t[i+1]) * bspline(t, i+1, k-1)(x)

        g(x) = f1(x) + f2(x)
        return g
    end
end

"""
    Boundary spline C3()
    Used in add_control_pts()
"""
function C3()
    f(x) = (1-x)^3 * bspline(0:1:1, 1, 0)(x)
    return f
end

"""
    Boundary spline C2()
    Used in add_control_pts()
"""
function C2()
    f(x) = x*(1-x)*(2-3x/2)*bspline(0:1:1, 1, 0)(x) + (2-x)^2*bspline(0:1:2, 1, 1)(x)/4
    return f
end

"""
    Boundary spline C1()
    Used in add_control_pts()
"""
function C1()
    f(x) = (x^2)*(1-x)*bspline(0:1:1, 1, 0)(x)/2 + x*(2-x)*bspline(0:1:2, 1, 1)(x)/4 + (3-x)*bspline(0:1:3, 1, 2)(x)/3
    return f
end


"""
    tester function
"""
function test_boundary_splines()
    x = [0:0.001:1]
    plot([0:0.0001:1], C3())
    plot!([0:0.0001:2], C2())
    plot!([0:0.0001:3], C1())
end


"""
    Creates a fitted B-spline around a set of "pts" (1d array), of the B-spline "B",
        pts of length "N", B-spline of degree "k"-1, of iteration round "iter"

    Notes:  As this is a recursive function, the parameter "iter" is optional and relays
        to the function itself
            Uses boundary splines (C1, C2, C3)
"""
function add_control_pts(pts, B, N, k, iter=N)
    if iter == 0
        f(x) = 0
        return f
    else
        n = N - 3
        g(x) = 
            if iter == 1
                pts[iter]*C3()(n*x)
            elseif iter == 2
                pts[iter]*C2()(n*x)
            elseif iter == 3   
                pts[iter]*C1()(n*x)

            elseif iter == N-2
                pts[iter]*C1()(n-n*x)
            elseif iter == N-1
                pts[iter]*C2()(n-n*x)
            elseif iter == N
                pts[iter]*C3()(n-n*x)
            else
                pts[iter]*B(n*x - (iter-3) + 1)
            end

                                    # decrease the iteration
        f = add_control_pts(pts, B, N, k, iter-1)     
        h(x) = f(x) + g(x)
        return h
    end
end