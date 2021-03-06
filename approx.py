import numpy as np
import matplotlib.pyplot as plt

def relaxation(f, g, exact_root, first_point, iterations=100, stop_condition=-1, args=[]):
    xr = np.zeros(iterations + 1)
    xr[0] = first_point
    error = np.full(iterations, np.inf)
    
    i = 0
    while (i < iterations) and (min(error) > stop_condition):
        xr[i+1] = g(xr[i])
        error[i] = abs(xr[i] - exact_root)
        
        points = xr[1:]
        
        i += 1
    
    return points[:i-1], error[:i-1], exact_root

def newton(f, exact_root, first_point, second_point, iterations=100, stop_condition=-1, args=[]):  
    xr = np.zeros(iterations + 2)
    derivative = np.zeros(iterations + 1)
    error = np.full(iterations, np.inf)
    
    xr[0], xr[1] = first_point, second_point
    first_derivative = (f(xr[1], *args) - f(xr[0], *args)) / (xr[1] - xr[0])
    derivative[0] = first_derivative
    
    i = 0
    while (i < iterations) and (min(error) > stop_condition):
        xr[i+2] = xr[i+1] - f(xr[i+1], *args)/derivative[i]
        derivative[i+1] = (f(xr[i+1], *args) - f(xr[i], *args)) / (xr[i+1] - xr[i])
        
        if np.isnan(xr[i+2]) == True:
            break
        
        error[i] = abs(xr[i] - exact_root)
            
        points = xr[2:]
        
        i += 1
        
    return points[:i], error[:i], exact_root

def secant(f, exact_root, first_point, second_point, iterations=100, stop_condition=-1, args=[]):
    xr = np.zeros(iterations+2)
    xr[0], xr[1] = first_point, second_point
    error = np.full(iterations, np.inf)
    
    i = 0
    while (i < iterations) and (min(error) > stop_condition):
        xr[i+2] = xr[i+1] - f(xr[i+1], *args)*(xr[i+1]-xr[i])/(f(xr[i+1], *args)-f(xr[i], *args))

        error[i] = abs(xr[i+2] - exact_root)

        points = xr[2:]
        
        i += 1
        
    return points[:i], error[:i], exact_root

def bolzano(f, exact_root, first_point, second_point, iterations=100, stop_condition=-1, args=[]):
    middle_point = (first_point + second_point)/2
    a = first_point
    b = second_point  
    c = np.zeros(iterations + 1)
    c[0] = middle_point

    sigma = (b-a)/2
    error = np.full(iterations, np.inf)
    
    i = 0
    while (i < iterations) and (min(error) > stop_condition):

        if f(a, *args)*f(c[i], *args) < 0:
            sigma = (b-a)/2
            c[i+1] = (a + c[i])/2
            a = c[i] - sigma
            b = c[i]

        elif f(c[i], *args)*f(b, *args) < 0:
            sigma = (b-a)/2
            c[i+1] = (c[i] + b)/2
            b = c[i] + sigma
            a = c[i]

        else:
            return points[:i], error[:i], exact_root

        error[i] = abs(c[i+1] - exact_root)
        
        points = c[1:]
        
        i += 1
        
    return points[:i-1], error[:i-1], exact_root

def false_position(f, exact_root, first_point, second_point, iterations=100, stop_condition=-1, args=[]):
    a = first_point
    b = second_point
    middle_point = (a*f(b, *args)-b*f(a, *args))/(f(b, *args)-f(a, *args))
    c = np.zeros(iterations + 1)
    c[0] = middle_point
    error = np.full(iterations, np.inf)
    
    i = 0
    while (i < iterations) and (min(error) > stop_condition):

        if f(a, *args)*f(c[i], *args) < 0:
            c[i+1] = (a*f(c[i], *args)-c[i]*f(a, *args))/(f(c[i], *args)-f(a, *args))
            b = c[i]

        elif f(c[i], *args)*f(b, *args) < 0:
            c[i+1] = (c[i]*f(b, *args)-b*f(c[i], *args))/(f(b, *args)-f(c[i], *args))
            a = c[i]

        else:
            return points, error, exact_root

        error[i] = abs(c[i+1] - exact_root)

        points = c[1:]
        
        i += 1
    
    return points[:i-1], error[:i-1], exact_root

def plot(title, points, f, exact_root, args=[]):
    x = np.linspace(exact_root-1, exact_root + 1, 100)
    plt.title(title)
    plt.plot([exact_root], [f(exact_root, *args)], color='k', marker='+', markersize=20, label='Root')
    plt.plot(x, f(x, *args), label='$f(x)$')
    plt.plot(points, f(points, *args), marker='o', linestyle='dashed', label='Iterations')
    plt.xlabel('$x$')
    plt.ylabel('$f(x)$')
    plt.legend(loc='best')
    
    plt.grid(linestyle='dashed')
    
    return None

def error_plot(title, error):
    iterations = len(error)
    
    x = np.linspace(1,iterations,iterations)
    
    plt.title('Error plot of '+title)
    plt.plot(x, error, label='Error')
    plt.yscale('log')
    plt.xlabel('Iterations')
    plt.ylabel('Error')
    plt.legend(loc='best')
    plt.grid(linestyle='dashed')
    
    return None