using Plots

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
    Creates a fitted B-spline around a set of "pts" (1d array), of the B-spline "B",
        pts of length "N", B-spline of degree "k"-1, of iteration round "iter"

    Note: as this is a recursive function, the parameter "iter" is optional and relys
        to the function itself
"""
function add_control_pts(pts, B, N, k, iter=N)
    if iter == 0 #-1
        f(x) = 0
        return f
    else
                                          # decrease the iteration
        f = add_control_pts(pts, B, N, k, iter-1)     
        g(x) = pts[iter]*B(N*x - iter + 2)    # pts[iter + 1]  due to julia indexing starting at 1    
        h(x) = f(x) + g(x)
        return h
    end

end


"""
    Draws a curve of degree "k"-1, with "N" points located in "x_pts", "y_pts"
"""
function run(k, N, x_pts, y_pts)
    X = Matrix{Float64}(undef, N, N)
    Y = Matrix{Float64}(undef, N, N)

    X .= x_pts
    Y .= y_pts

    knots = 0 : 1 : k+1
    B = bspline(knots, 1, k)
    
    Bx = add_control_pts(X, B, N, k)
    By = add_control_pts(Y, B, N, k)

    x = 0:0.0001:1

    scatter(X, Y)
    plot!(Bx.(x), By.(x), legend=false)
end

function a()
    run(3, 9, [7.2, 3.35, -7.05, -11.1, 2.05, 9.2, 7.25, 6.95, 15.25], [7,   10,    6.4, -8.6,   -6.5, 9.8, 3.7, -13.1, -0.8])

end

function example1()
    run(3, 8, [0, 1, 2, 3, 4, 5, 6, 5], [.2, 0, .2, 1, 1.3, 1.3, 1, .8])

end

# example1() 
# a()     #comment to view example1
