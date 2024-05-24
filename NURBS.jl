using Plots

function bspline(t, i, k)
    if k == 1
        f(x) = t[i] <= x < t[i+1] ? 1 : 0
        return f
    else
        f1(x) = (x-t[i])/(t[i+k-1]-t[i]) * bspline(t, i, k-1)(x)
        f2(x) = (t[i+k]-x)/(t[i+k]-t[i+1]) * bspline(t, i+1, k-1)(x)

        g(x) = f1(x) + f2(x)
        return g
    end
end
# function bspline(t, i, k)
#     if k == 0
#         f(x) = t[i] <= x < t[i+1] ? 1 : 0
#         return f
#     else
#         f1(x) = (x-t[i])/(t[i+k]-t[i]) * bspline(t, i, k-1)(x)
#         f2(x) = (t[i+k+1]-x)/(t[i+k+1]-t[i+1]) * bspline(t, i+1, k-1)(x)

#         g(x) = f1(x) + f2(x)
    
#         return g
#     end
# end


""" degree n, k number of knots
"""
function numerator(knots, k, w, pts, iter=length(w))
    if iter == 0
        a(x) = 0
        return a
    else
        f(x) = w[iter]*pts[iter]*bspline(knots, iter, k)(x)
        g = numerator(knots, k, w, pts, iter-1)
        h(x) = f(x) + g(x)
        return h
    end
end

function denominator(knots, k, w, iter=length(w))
    if iter == 0
        a(x) = 0
        return a
    else
        f(x) = w[iter]*bspline(knots, iter, k)(x)
        g = denominator(knots, k, w, iter-1)
        h(x) = f(x) + g(x)
        return h
    end
end
    

function Q(knots, k, w, xpts, ypts)
    N = length(w)
    b = denominator(knots, k, w)

    ax = numerator(knots, k, w, xpts)
    Qx(x) = ax(x)/b(x) 

    ay = numerator(knots, k, w, ypts)
    Qy(x) = ay(x)/b(x)
    
    return (Qx, Qy)
end

function graph()
    #knots = [0, 0, 0, pi/2, pi/2, pi, pi, 3pi/2, 3pi/2, 2pi, 2pi, 2pi]
    #knots = [0.0, 0.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    knots = [0.0, 1.0, 2.0, 3.0, 3.0, 5.0, 6.0, 7.0]


    b = bspline(knots, 1, 3)
    val = 0:0.01:2pi
    p = plot(val, b.(val))
    # display(p)
end


"""page 63"""
function test()
    # knots = [0, 0, 0, 1/4, 1/4, 1/2, 1/2, 3/4, 3/4, 1, 1, 1]
    # x = [1, sqrt(2)/2, 0, -sqrt(2)/2, -1, -sqrt(2)/2, 0, sqrt(2)/2, 1] 
    # y = [0, sqrt(2)/2, 1, sqrt(2)/2, 0, -sqrt(2)/2, -1, -sqrt(2)/2, 0]
    # w = [1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1]
#    knots = [0, 0, 0, pi/2, pi/2, pi, pi, 3pi/2, 3pi/2, 2pi, 2pi, 2pi]
    knots = [0, pi/8, pi/4, 3pi/8, pi/2, 5pi/8, 3pi/4, 7pi/8, pi, 9pi/8, 5pi/4, 11pi/8, 3pi/2, 13pi/8, 7pi/4, 15pi/8, 2pi]
    x = [1, 1, 0, -1, -1, -1, 0, 1, 1] 
    y = [0, 1, 1, 1, 0, -1, -1, -1, 0]
    w = [1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1]

    k = 7

    Qx, Qy = Q(knots, k, w, x, y)

    val = 0:0.01:2pi
    p = plot(Qx.(val), Qy.(val))
    # scatter!(x, y)

    # println(Qx(0), Qy(0))

    # display(p)
end

test()
#graph()