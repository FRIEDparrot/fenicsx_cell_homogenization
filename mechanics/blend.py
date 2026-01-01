import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import logm, expm

def log_euclidean_blend(A, B, steps=100):
    """
    Log-Euclidean blend of two matrices A and B.
    Args:
        A:
        B:
        steps: number of steps to perform the blend. Default is 100.
    Returns:
    """
    t_spec = np.linspace(0, 1, steps)
    def s(t):
        # cubic easing function
        return  3 * t**2 - 2 * t**3
    result = np.array([expm(logm(A) * (1 - s(t)) + logm(B) * s(t)) for t in t_spec])
    return result  # log-euclidean blend


if __name__ == "__main__":
    a = np.array([[100]])
    b = np.array([[300]])
    steps = 100
    x = np.linspace(0.01, 1, steps)
    y = log_euclidean_blend(a, b, steps)
    plt.plot(x, y.reshape(-1))
    plt.show()
