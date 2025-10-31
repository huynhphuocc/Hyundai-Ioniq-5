# GA Crash Test Optimization - Mô hình 2DOF

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import random

# ----- 1. Mô hình mô phỏng crash 2DOF -----
def simulate_crash(params, t_span, v0):
    k1, c1, k2, c2 = params
    m1, m2 = 1500, 75  # Khối lượng xe và người (kg)

    def model(t, y):
        x1, v1, x2, v2 = y
        dx1dt = v1
        dx2dt = v2
        dv1dt = (-k1 * x1 - c1 * v1 + k2 * (x2 - x1) + c2 * (v2 - v1)) / m1
        dv2dt = (-k2 * (x2 - x1) - c2 * (v2 - v1)) / m2
        return [dx1dt, dv1dt, dx2dt, dv2dt]

    y0 = [0, -v0, 0, -v0]  # Cùng tốc độ ban đầu tiến về phía trước
    sol = solve_ivp(model, [t_span[0], t_span[-1]], y0, t_eval=t_span, method='RK45')
    return sol.y  # Trả về [x1, v1, x2, v2]

# ----- 2. Tạo dữ liệu ground truth giả lập -----
def generate_ground_truth():
    t = np.linspace(0, 0.2, 1000)
    v0 = 13.89  # 50 km/h
    true_params = [22000, 1500, 9000, 1200]  # k1, c1, k2, c2
    y = simulate_crash(true_params, t, v0)
    return t, y, true_params

# ----- 3. Hàm fitness tính sai số -----
def fitness(params, t, y_true, v0):
    try:
        y_pred = simulate_crash(params, t, v0)
        error = np.mean((y_pred[0] - y_true[0])**2 + (y_pred[2] - y_true[2])**2)  # so sánh x1, x2
        return error
    except:
        return 1e9  # Phạt nếu thông số gây mô hình nổ

# ----- 4. Thuật toán GA đơn giản -----
def run_ga(t, y_true, v0, generations=50, pop_size=20):
    bounds = [(1e3, 3e4), (500, 3000), (1e3, 1e4), (300, 3000)]
    population = [np.array([random.uniform(*b) for b in bounds]) for _ in range(pop_size)]

    for gen in range(generations):
        scored = [(ind, fitness(ind, t, y_true, v0)) for ind in population]
        scored.sort(key=lambda x: x[1])
        elites = [x[0] for x in scored[:pop_size//2]]

        # Tạo thế hệ mới bằng lai ghép + đột biến
        children = []
        while len(children) < pop_size:
            p1, p2 = random.sample(elites, 2)
            child = [(a + b) / 2 for a, b in zip(p1, p2)]
            child = [c * random.uniform(0.9, 1.1) for c in child]  # mutation
            child = [np.clip(c, *bounds[i]) for i, c in enumerate(child)]
            children.append(np.array(child))

        population = children
        print(f"Gen {gen}: Best Fitness = {scored[0][1]:.4f}")

    return scored[0][0]  # Trả về bộ tham số tốt nhất

# ----- 5. Chạy toàn bộ pipeline -----
t, y_true, true_params = generate_ground_truth()
v0 = 13.89
best_params = run_ga(t, y_true, v0)
print("\nTham số đúng:", true_params)
print("Tham số tìm được:", best_params)

# ----- 6. Vẽ kết quả -----
y_sim = simulate_crash(best_params, t, v0)
plt.plot(t, y_true[0], label='x1 True', linestyle='--')
plt.plot(t, y_sim[0], label='x1 GA')
plt.plot(t, y_true[2], label='x2 True', linestyle='--')
plt.plot(t, y_sim[2], label='x2 GA')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend()
plt.title('So sánh mô phỏng và dữ liệu chuẩn')
plt.grid(True)
plt.show()
