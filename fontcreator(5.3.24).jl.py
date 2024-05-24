import matplotlib.pyplot as plt
import numpy as np

def graph(points):
    # start, stop, step
    x = np.arange(0., 5., 0.2)
    
    plt.plot(x, np.piecewise(x, [x  == 2, x != 2],  # when x
                                [x/x, 2]))      # func 
    plt.show()


def add_splines(k, funclist1, funclist2, condlist1, condlist2):
    newfunclist = [funclist1[0]]
    
    for i in range(k):
        newfunclist.append(np.add(funclist1[i+1], funclist2[i]))
    newfunclist.append(funclist2[-1])

    condlist1.append(condlist2[-1])
    newcondlist = condlist1

    return newfunclist, newcondlist

def make_spline(x, t, i, k) :
    """
        Returns "funcList", "condList"
    """

    if k == 0:
        return [1, 0], [(t[i] <= x)*(x < t[i + 1])]
    else:
        func1, cond1 = make_spline(x, t, i, k-1)
        func2, cond2 = make_spline(x, t, i+1, k-1)
    

        for j in range(k):
             func1[j] = func1[j] #* (x-t[i])/(t[i+k]-t[i])
        # for j in range(k):
        #     func2[j] = np.multiply(func2[j], (t[i+k+1]-x)/(t[i+k+1]-t[i-1]))


        return add_splines(k, func1, func2, cond1, cond2)

if __name__ == '__main__':
    #x = np.array([7.2, 3.35, -7.05, -11.1, 2.05, 9.2, 7.25, 6.95, 15.25])
    #y = np.array([7,   10,    6.4, -8.6,   -6.5, 9.8, 3.7, -13.1, -0.8])
    
    x = np.arange(0., 5., .1)

    #print([(1 <= x)*(x < 2), x < 0], "end")

    # funclist, condlist = make_spline(x, [1, 2, 3, 4], 0, 1)

    # newfunclist = []
    
    # plt.plot(x, np.piecewise(x, condlist, funclist)) 
    # plt.show()
    graph(x)

