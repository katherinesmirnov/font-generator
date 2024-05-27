using Plots

"""
References: 
    (1) Apple Tech Journal: https://vintageapple.org/develop/pdf/develop-25_9603_March_1996.pdf#page=50
    (2) Wiki: https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline#cite_note-nurbs-book-11
    (3) Piegl, Tiller: https://argos.vu/wp-content/uploads/2019/03/Springer-TheNURBSBook.pdf


    Other help:
        https://computergraphics.stackexchange.com/questions/7908/drawing-nurbs-knots-multiplicity-divide-by-zero
        http://web.cs.wpi.edu/~matt/courses/cs563/talks/nurbs.html

"""

function bspline(t, i, k)
    if k == 0
        f(x) = t[i] <= x < t[i+1] ? 1 : 0
        return f
    else
        a(x) = (t[i+k]-t[i]) == 0 ? 0 : (x-t[i])/(t[i+k]-t[i])
        b(x) = (t[i+k+1]-t[i+1]) == 0 ? 0 : (t[i+k+1]-x)/(t[i+k+1]-t[i+1])

        f1(x) =  a(x) * bspline(t, i, k-1)(x) 
        f2(x) =  b(x) * bspline(t, i+1, k-1)(x) 

        g(x) = f1(x) + f2(x)

        return g
    end
end


function Q(knots, k, w, xpts, ypts)
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

    b = denominator(knots, k, w)

    ax = numerator(knots, k, w, xpts)
    Qx(x) = ax(x)/b(x) 

    ay = numerator(knots, k, w, ypts)
    Qy(x) = ay(x)/b(x)

    return (Qx, Qy)
end


"""
    knots u, v
    deg "p" in u direction
    def "q" in v direction
    w: 2d array of weights (n x m)
"""
function Qsurface1(u_knots, v_knots, u_degree, v_degree, weights, coords)
    n = length(weights)
    m = length(weights[1])

    println(n, " x ", m)

    function R_numerator(i, j)
        f1(u, v) = bspline(u_knots, i, u_degree)(u)
        f2(u, v) = bspline(v_knots, j, v_degree)(v)
        f(u, v) = f1(u, v) * f2(u, v) * weights[i][j]
        return f
    end

    function R_denominator(i, j, k=n)
        function inner(l=m)
            if l == 0
                g(u, v) = 0
                return f
            else
                f1(u, v) = bspline(v_knots, l, v_degree)(v)*weights[k][l]
                f2(u, v) = inner(l -1)(u, v)
                f(u, v) = f1(u, v) + f2(u, v)
                return f

                # f(u, v) = inner(l-1)(u, v)
                # return f
                
            end
        end

        if k == 0
            g(u, v) = 0
            return g
        else
            f1(u, v) = bspline(u_knots, k, u_degree)(u)
            f2(u, v) = inner()(u, v)

            f3(u, v) = R_denominator(i, j, k-1)(u, v)
            f(u, v) = f1(u, v)*f2(u, v) + f3(u, v)
            return f
        end        
    end

    function R(i, j)
        f1(u, v) = R_numerator(i, j)(u, v)
        f2(u, v) = R_denominator(i, j)(u, v)
        f(u, v) = f1(u, v)/f2(u, v)
        return f
    end

    function S(i = m)
        function inner(j=m)
            if j == 0
                g(u, v) = 0
                return g
            else
                f1(u, v) = R(i, j)(u, v)*coords[i][j]
                f2(u, v) = inner(j-1)(u, v)
                f(u, v) = f1(u, v) + f2(u, v)
                return f
                # f(u, v) = inner(j-1)(u, v)
                # return f
            end
        end

        if i == 0
            g(u, v) = 0
            return g
        else
            f1(u, v) = inner()(u, v)
            f2(u, v) = S(i - 1)(u, v)
            f(u, v) = f1(u, v) + f2(u, v)
            return f
            # f(u, v) = S(i - 1)(u, v)
            # return f
        end
    end

    f(u, v) = S()(u, v) 
    return f
end

function Qsurface(u, v, u_knots, v_knots, u_degree, v_degree, weights, coords)
    n = length(weights)
    m = length(weights[1])
    
    function R(i, j)
        numerator = bspline(u_knots, i, u_degree)(u) * bspline(v_knots, j, v_degree)(v) * weights[i][j]
        if numerator == 0.0
            return 0.0
        end

        denominator = 0.0
        for k in 1:1:n
            for l in 1:1:m
                denominator += bspline(u_knots, k, u_degree)(u) * bspline(v_knots, l, v_degree)(v) * weights[k][l]
            end
        end

        return numerator/denominator
    end

    
    function S()
        f = [0.0, 0.0, 0.0]
        for i in 1:1:n
            for j in 1:1:m
                f += R(i, j)*coords[i][j]
            end
        end 
        return f
    end

    
    return S()
end

""" Figure 12 from (1)
"""
function graph12()
    knots = [0.0, 0.0, 0.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    b = bspline(knots, 1, 2)
    val = 0:0.01:7
    p = plot(val, b.(val))

    plot!(val, bspline(knots, 2, 2).(val))
    plot!(val, bspline(knots, 3, 2).(val))
    plot!(val, bspline(knots, 4, 2).(val))
    plot!(val, bspline(knots, 5, 2).(val))

    display(p)
end

""" Figure 14 from (1)
"""
function graph14()
    knots = [.0, 1.0, 2.0, 3.0, 3.0, 5.0, 6.0, 7.0]

    b = bspline(knots, 1, 2)
    val = 0:0.01:7
    p = plot(val, b.(val))

    plot!(val, bspline(knots, 2, 2).(val))
    plot!(val, bspline(knots, 3, 2).(val))
    plot!(val, bspline(knots, 4, 2).(val))
    plot!(val, bspline(knots, 5, 2).(val))

    display(p)
end


""" Knots from (2) Wiki """
function circle()
    knots = [0, 0, 0, pi/2, pi/2, pi, pi, 3pi/2, 3pi/2, 2pi, 2pi, 2pi]
    x = [1, 1, 0, -1, -1, -1, 0, 1, 1] 
    y = [0, 1, 1, 1, 0, -1, -1, -1, 0]
    w = [1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1, sqrt(2)/2, 1]

    k = length(knots) - length(x) - 1
    println("Degree: ", k+1)

    Qx, Qy = Q(knots, k, w, x, y)

    val = 0:0.001:2pi
    p = plot(Qx.(val), Qy.(val))
    scatter!(x, y)

    display(p)
end


function sphere()
    u_deg = 2
    v_deg = 2 

    u_knots = Float64[0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5]
    v_knots = Float64[0, 0, 0, 1, 2, 3, 3, 3]

    u_knots ./= 5.0
    v_knots ./= 3.0

    control_points = [
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 2.0, 4.0], [0.0, 3.0, 2.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [2.0, 3.0, 4.0], [2.0, 4.0, 2.0], [2.0, 3.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [4.0, 2.0, 4.0], [4.0, 3.0, 2.0], [4.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    ]

    weights = [
        [0.0 0.0 0.0 0.0 0.0],
        [0.0 0.0 0.0 0.0 0.0],
        [0.0 1.0 2.0 1.0 0.0],
        [0.0 2.0 6.0 2.0 0.0],
        [0.0 1.0 2.0 1.0 0.0],
        [0.0 0.0 0.0 0.0 0.0],
        [0.0 0.0 0.0 0.0 0.0],
        [0.0 0.0 0.0 0.0 0.0]
    ]

    t = range(0, 1, length=10)

    @assert length(weights) == length(control_points)
    @assert length(weights[1]) == length(control_points[1])


    # S = Qsurface1(u_knots, v_knots, u_deg, v_deg, weights, control_points)
    # p = plot(S.(t, t), S.(t, t), st =:surface, legend=false)

    # x_vals = []
    # y_vals = []
    # z_vals = []

    # for u in t
    #     ux_vals = Float64[]
    #     uy_vals = Float64[]
    #     uz_vals = Float64[]

    #     for v in t
    #         x, y, z = Qsurface(u, v, u_knots, v_knots, u_deg, v_deg, weights, control_points)
    #         append!(ux_vals, x)
    #         append!(uy_vals, y)
    #         append!(uz_vals, z)

    #     end
    #     push!(x_vals, ux_vals)
    #     push!(y_vals, uy_vals)
    #     push!(z_vals, uz_vals)
    # end

    surface_pts = [Qsurface(u, v, u_knots, v_knots, u_deg, v_deg, weights, control_points) for u in t, v in t]

    x_vals = [lst[1] for lst in surface_pts]
    y_vals = [lst[2] for lst in surface_pts]
    z_vals = [lst[3] for lst in surface_pts]

    print(x_vals)
    #p = plot(x_vals[1], z_vals[1], st =:surface, legend=false)
    p = plot3d(t, t, z_vals, st =:surface, legend=false)

    
    display(p)
end

sphere()