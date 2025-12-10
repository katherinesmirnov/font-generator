include("fontcreator.jl")
using Roots

small_letter = 9
tall_letter = 15
connecting_height = 5

lettersize = Dict(
    'a' => small_letter,
    'b' => tall_letter,
    'c' => small_letter,
    'd' => tall_letter,
    'e' => small_letter,
    'f' => tall_letter,
    'g' => tall_letter,
    'h' => tall_letter,
    'i' => small_letter,
    'j' => tall_letter,
    'k' => tall_letter,
    'l' => tall_letter,
    'm' => small_letter,
    'n' => small_letter,
    'o' => small_letter,
    'p' => tall_letter,
    'q' => tall_letter,
    'r' => small_letter,
    's' => small_letter,
    't' => tall_letter,
    'u' => small_letter,
    'v' => small_letter,
    'w' => small_letter,
    'x' => small_letter,
    'y' => tall_letter,
    'z' => tall_letter,
)

offsetY = Dict(
    'b' => -1,
    'f' => -2,
    'g' => -8,
    'j' => -6,
    'p' => -6,
    'q' => -6,
    'y' => -7,
    'z' => -7,
)


#https://public.vrac.iastate.edu/~oliver/courses/me625/week5b.pdf

function getExtrema(Bx, By)
    maxI = 0 
    maxY = -1000

    minI = 0 
    minY = 1000
    for i in 0:0.001:1
        if (By(i) > maxY)
            maxY = By(i)
            maxI = i 
        elseif (By(i) < minY)
            minY = By(i)
            minI = i 
        end
    end
    return minI, maxI
end

function shiftValues(yPts, minI)
    for i in 1:1:length(yPts)
        yPts[i] -= minI
    end
end

function scaleValues(xPts, yPts, height, desiredHeight)
    for i in 1:1:length(xPts)
        xPts[i] *= (desiredHeight)/(height)
        yPts[i] *= (desiredHeight)/(height)

    end
end

function addOffset(yPts, offset)
    for i in 1:1:length(yPts)
        yPts[i] += offset
    end
end

function modifyEndponts(yPts)
    yPts[1] = connecting_height
    yPts[end] = connecting_height
end

function scale(letter, (x, y))
    Bx = add_control_pts(x, bspline( 0 :1 : 3+1, 1, 3), length(x), 3)
    By = add_control_pts(y, bspline( 0 :1 : 3+1, 1, 3), length(y), 3)
    i0, i1 = getExtrema(Bx, By)

    shiftValues(y, By(i0))
    scaleValues(x, y, By(i1)-By(i0), lettersize[letter])
    if haskey(offsetY, letter)
        addOffset(y, offsetY[letter]) 
    end
    modifyEndponts(y)

    Bx = add_control_pts(x, bspline( 0 :1 : 3+1, 1, 3), length(x), 3)
end    

function testalpha()
    str = ""
    for key in sort(collect(keys(coordinates)))
        str *= key
        scale(key, coordinates[key])
    end

    writeSentence(str)    
end


function calcBprime(pts, k, N, iter = N-1)
    if iter == 0 #-1
        h(x) = 0
        return h
    else 
        K = k+1
        # if iter == 0
        #     ((K)/(iter+K+1 - iter+1))*(pts[iter+1])
        # else
        #     ((K)/(iter+K+1 - iter+1))*(pts[iter+1] - pts[iter])
        # end

        
        # knots = iter+1 : 1: iter+k+2
        # B = bspline(knots, 1, k)   
        knots = 0 : 1 : 4       
        knots = vcat([0, 0], 0 : 1/(N) : 1)
        knots = vcat(knots[k], [1, 1])
        println(knots)
        B = bspline(knots, iter, k)   

        # knots = 0 : 1/(N+k-1) : (N+k)/(N+k-1)
        # B = bspline(knots, iter, k)   

        a(x) = B(x)*((K)/(knots[iter+K+1] - knots[iter+1]))*(pts[iter+1] - pts[iter])
        b = calcBprime(pts, k, N, iter -1)
        c(x) = a(x) + b(x) 
        return c
    end
end

function test()
    x, y = coordinates['b']
    #x, y = [-22.5,-9.7,-2.9, 0, 3, 4], [ 25.6, 16.6, 13.3, 10, 20, 30]
    k = 2
    t = 0:0.01:1

    Bx = add_control_pts(x, bspline( 0 :1 : 3+1, 1, 3), length(x), 3)
    By = add_control_pts(y, bspline( 0 :1 : 3+1, 1, 3), length(y), 3)
    #scatter!(x, y)

    
    Bxprime = calcBprime(x, k, length(x))
    Byprime = calcBprime(y, k, length(x))
    
    Bprime(x) = Byprime(x)/Bxprime(x)
  
    #plot(t, Bxprime.(t), legend=false, ylim=(-10, 10))
    #plot!(t, Byprime.(t), legend=false)
  
    

    #p = plot(Bx.(t), By.(t), legend=false) #, ylim=(-10, 10))
    p = plot(t, Bprime.(t), legend=false) #, ylim=(-10, 10))  


    # z = find_zeros(Bprime, 0, 1) 
    # max = 10000
    # minsize = -10000
    # for zero in z
    #     #scatter!([Bx(zero)], [By(zero)])
    # end
    
    display(p)
end

function printDict()
    for key in sort(collect(keys(coordinates)))
        print("'", key, "'", " => ")
        println("(", coordinates[key][1], ",")
        println(coordinates[key][2], "),")
    end
end

function testExtrema()
    x, y = coordinates['b']
    Bx = add_control_pts(x, bspline( 0 :1 : 3+1, 1, 3), length(x), 3)
    By = add_control_pts(y, bspline( 0 :1 : 3+1, 1, 3), length(y), 3)
    i0, i1 = getExtrema(Bx, By)

    t = 0:0.01:1
    p = plot(Bx.(t), By.(t), legend=false)
    scatter!([Bx(i0)], [By(i0)])
    scatter!([Bx(i1)], [By(i1)])
    display(p)
end

#testalpha()
#write("abc")
#printDict()
testExtrema()