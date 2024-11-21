import numpy as np
from scipy.optimize import minimize, differential_evolution
from functools import lru_cache
import concurrent.futures

# Example nonlinear function (replace this with your actual slow function)
@lru_cache(maxsize=None)
def nonlinear_function(x):
    x1, x2, x3, x4, x5, x6 = x
    return (np.sin(x1) + np.cos(x2) + np.tan(x3) +
            np.exp(x4) - np.log(np.abs(x5) + 1) + x6**2)

# Parallel wrapper function
def parallel_function(x):
    return nonlinear_function(tuple(x))  # Convert to tuple for caching

# Global optimization step using Differential Evolution (good for expensive functions)
bounds = [(-10, 10)] * 6  # Adjust bounds as needed

def optimize_function():
    # Step 1: Use Differential Evolution to find a good initial point
    print("Running global optimization with Differential Evolution...")
    global_result = differential_evolution(parallel_function, bounds, maxiter=10, workers=-1)
    
    print(f"Global optimization result: {global_result.x}")
    
    # Step 2: Use a local optimizer to refine the result
    print("Refining solution with local optimization...")
    local_result = minimize(parallel_function, global_result.x, method='Nelder-Mead')
    
    return local_result

if __name__ == "__main__":
    result = optimize_function()
    
    # Display the results
    print("\nOptimization Results:")
    print(f"Optimal input variables: {result.x}")
    print(f"Minimum function value: {result.fun}")
    print(f"Number of iterations: {result.nit}")
    print(f"Success status: {result.success}")
    print(f"Message: {result.message}")