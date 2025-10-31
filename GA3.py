# GA Crash Test Optimization - 2DOF Model

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random

# ----- 1. 2DOF Crash Simulation Model -----
def simulate_crash(params, t_span, v0):
    k1, c1, k2, c2 = params
    m1, m2 = 1500, 75  # Mass of vehicle and occupant (kg)

    def model(t, y):
        x1, v1, x2, v2 = y
        dx1dt = v1
        dx2dt = v2
        dv1dt = (-k1 * x1 - c1 * v1 + k2 * (x2 - x1) + c2 * (v2 - v1)) / m1
        dv2dt = (-k2 * (x2 - x1) - c2 * (v2 - v1)) / m2
        return [dx1dt, dv1dt, dx2dt, dv2dt]

    y0 = [0, -v0, 0, -v0]  # Initial velocity toward impact
    sol = solve_ivp(model, [t_span[0], t_span[-1]], y0, t_eval=t_span, method='RK45')
    return sol.y  # Returns [x1, v1, x2, v2]

# ----- 2. Generate Simulated Ground Truth Data -----
def generate_ground_truth():
    t = np.linspace(0, 0.2, 1000)
    v0 = 13.89  # 50 km/h
    true_params = [22000, 1500, 9000, 1200]  # k1, c1, k2, c2
    y = simulate_crash(true_params, t, v0)
    return t, y, true_params

# ----- 3. Fitness Function for Error Calculation -----
def fitness(params, t, y_true, v0):
    try:
        y_pred = simulate_crash(params, t, v0)
        error = np.mean((y_pred[0] - y_true[0])**2 + (y_pred[2] - y_true[2])**2)  # Compare x1 and x2
        return error
    except:
        return 1e9  # Penalize for unstable simulations

# ----- 4. Basic Genetic Algorithm -----
def run_ga(t, y_true, v0, generations=50, pop_size=20):
    bounds = [(1e3, 3e4), (500, 3000), (1e3, 1e4), (300, 3000)]
    population = [np.array([random.uniform(*b) for b in bounds]) for _ in range(pop_size)]

    for gen in range(generations):
        scored = [(ind, fitness(ind, t, y_true, v0)) for ind in population]
        scored.sort(key=lambda x: x[1])
        elites = [x[0] for x in scored[:pop_size//2]]

        # Create new generation via crossover and mutation
        children = []
        while len(children) < pop_size:
            p1, p2 = random.sample(elites, 2)
            child = [(a + b) / 2 for a, b in zip(p1, p2)]
            child = [c * random.uniform(0.9, 1.1) for c in child]  # mutation
            child = [np.clip(c, *bounds[i]) for i, c in enumerate(child)]
            children.append(np.array(child))

        population = children
        print(f"Gen {gen}: Best Fitness = {scored[0][1]:.4f}")

    return scored[0][0]  # Return best parameter set

# ----- 5. Run Full Pipeline -----
t, y_true, true_params = generate_ground_truth()
v0 = 13.89
best_params = run_ga(t, y_true, v0)
print("\nGround truth parameters:", true_params)
print("Optimized parameters found by GA:", best_params)

# ----- 6. Plot Result -----
y_sim = simulate_crash(best_params, t, v0)
plt.plot(t, y_true[0], label='x1 Ground Truth', linestyle='--')
plt.plot(t, y_sim[0], label='x1 GA Output')
plt.plot(t, y_true[2], label='x2 Ground Truth', linestyle='--')
plt.plot(t, y_sim[2], label='x2 GA Output')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend()
plt.title('Comparison Between Simulation and Ground Truth')
plt.grid(True)
plt.show()
