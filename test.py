import numpy as np
import matplotlib.pyplot as plt

def basis_function(t, i, k, knots):
    """
    Compute the i-th B-spline basis function of degree k with knots at the given t.
    """
    if k == 0:
        return 1.0 if knots[i] <= t < knots[i + 1] else 0.0
    
    denom1 = knots[i + k] - knots[i]
    denom2 = knots[i + k + 1] - knots[i + 1]
    
    result = 0.0
    if denom1 != 0:
        result += (t - knots[i]) / denom1 * basis_function(t, i, k - 1, knots)
    if denom2 != 0:
        result += (knots[i + k + 1] - t) / denom2 * basis_function(t, i + 1, k - 1, knots)
    
    return result

def create_bspline(x, y, degree, num_points=1000):
    """
    Create a B-spline curve of the given degree using the given control points (x, y).
    """
    assert len(x) == len(y), "x and y must have the same length"
    assert len(x) >= degree + 1, "Not enough control points for the given degree"

    n = len(x) - 1  # number of control points
    m = n + degree + 1  # total number of knots
    knots = np.linspace(x[0], x[-1], m)  # uniform knot vector

    curve_x = np.linspace(x[0], x[-1], num_points)
    curve_y = np.zeros_like(curve_x)

    for i, t in enumerate(curve_x):
        for j in range(n):
            curve_y[i] += y[j] * basis_function(t, j, degree, knots)

    return curve_x, curve_y

def main():
    # Example data points

    x = np.array([7.2, 3.35, -7.05, -11.1, 2.05, 9.2, 7.25, 6.95, 15.25])
    y = np.array([7,   10,    6.4, -8.6,   -6.5, 9.8, 3.7, -13.1, -0.8])

    # Specify the degree of the B-spline
    degree = int(input("Enter the degree of the B-spline: "))

    # Create the B-spline curve
    curve_x, curve_y = create_bspline(x, y, degree)

    # Plot the original data points and the B-spline curve
    plt.figure()
    plt.plot(x, y, 'bo', label='Original points')
    plt.plot(curve_x, curve_y, 'r-', label=f'B-spline (degree {degree})')
    plt.legend()
    plt.title('B-spline interpolation (From Scratch)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
