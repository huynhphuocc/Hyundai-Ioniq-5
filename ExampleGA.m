%% ExampleGA.m  —  GA fit 2-DOF Vehicle–Occupant frontal crash
clear; clc; rng(0);

%% ===== 1) Inputs cơ bản =====
m_vehicle = 1905;                 % [kg] khối lượng xe
m1 = 0.65*m_vehicle;              % [kg] mass tương đương phần trước xe (0.6–0.75*m)
m2 = 75;                          % [kg] hành khách
v0 = 50/3.6;                      % [m/s] tốc độ ban đầu (50 km/h)

t_end = 0.18;                     % [s] thời lượng pulse ~ 180 ms
Nref  = 1801;                     % số mẫu tham chiếu (>=1000)
t_ref = linspace(0, t_end, Nref)';% lưới thời gian tham chiếu đều

%% ===== 2) Target pulse (a_meas): đổi tại đây nếu có dữ liệu thật =====
% DEMO: tạo pulse tam giác (thay bằng data thật là best)
a_peak_target_g = 35;             % [g] đỉnh mong muốn (ví dụ 35 g)
Tp = 0.09;                        % [s] độ dài pulse ~ 90 ms
a_meas = triangle_pulse(t_ref, a_peak_target_g*9.81, Tp);  % [m/s^2]

% Nếu có CSV thật (thời gian và gia tốc theo g), dùng:
% data = readmatrix('YourPulse.csv');  % cột1: t [s], cột2: a[g] hoặc a[m/s^2]
% t_raw = data(:,1); a_raw = data(:,2);
% a_meas = interp1(t_raw, a_raw*(abs(max(a_raw))<5)*9.81 + a_raw*(abs(max(a_raw))>=5), t_ref, 'linear', 'extrap');
% (dòng trên: nếu thấy trị g nhỏ, auto chuyển sang m/s^2; tuỳ file bạn chỉnh cho gọn)

%% ===== 3) GA bounds cho [k1 c1 k2 c2] =====
% k1,c1: phần xe (m1) — k2,c2: phần restraint (giữa m1–m2)
lb = [3e5,  2e4,  1e4,  1e3];
ub = [2e6,  1.5e5, 2e5, 5e4];
nvar = 4;

%% ===== 4) Objective function (RMSE gia tốc xe) =====
obj = @(p) fitness_2dof_align(p, t_ref, v0, m1, m2, a_meas);

%% ===== 5) GA options & run =====
if license('test','gads_toolbox') || exist('ga','file')
    opts = optimoptions('ga','PopulationSize',80,'MaxGenerations',120,...
        'Display','iter','UseParallel',false);
    [p_opt, fval] = ga(obj, nvar, [],[],[],[], lb, ub, [], opts);
else
    % Fallback: multi-start fminsearch nếu không có GA
    fprintf(' Không có Global Optimization Toolbox → dùng fminsearch multi-start.\n');
    best = inf; p_opt = [];
    for s = 1:12
        p0 = lb + (ub-lb).*rand(1,4);
        [p, f] = fminsearch(@(pp) fitness_2dof_align(pp, t_ref, v0, m1, m2, a_meas), ...
                            p0, optimset('Display','off','MaxFunEvals',8e3,'MaxIter',4e3));
        if f < best, best = f; p_opt = p; end
    end
    fval = best;
end
k1 = p_opt(1); c1 = p_opt(2); k2 = p_opt(3); c2 = p_opt(4);
fprintf('\n=== THÔNG SỐ TỐI ƯU ===\n');
fprintf('k1 = %.3e N/m,   c1 = %.3e N·s/m\n', k1, c1);
fprintf('k2 = %.3e N/m,   c2 = %.3e N·s/m\n', k2, c2);
fprintf('RMSE accel (SI) = %.3f m/s^2 (≈ %.2f g)\n', fval, fval/9.81);

%% ===== 6) Sim lại & plot =====
[t, a1, a2, x1, x2, v1, v2] = sim_2dof_interp(p_opt, t_ref, v0, m1, m2);

figure('Name','Vehicle–Occupant 2-DOF GA Fit','Color','w');
tiledlayout(3,1);

nexttile; 
plot(t_ref, a_meas/9.81, 'k--','LineWidth',1.2); hold on;
plot(t,     a1/9.81,      'b','LineWidth',1.5);
ylabel('a_{veh} [g]'); grid on; legend('target','sim'); title('Vehicle Accel');

nexttile;
plot(t, a2/9.81,'r','LineWidth',1.5);
ylabel('a_{occ} [g]'); grid on; title('Occupant Accel');

nexttile;
plot(t, x1,'b','LineWidth',1.3); hold on;
plot(t, x2,'r','LineWidth',1.3);
ylabel('x [m]'); xlabel('t [s]'); grid on; legend('x_1 (vehicle)','x_2 (occupant)');
title('Displacements');

%% ================= Local functions =================
function err = fitness_2dof_align(p, t_ref, v0, m1, m2, a_ref)
% Align chiều dài + penalty nếu ODE fail/NaN
    try
        [t, a1] = sim_vehicle_accel(p, t_ref, v0, m1, m2);
        % Interp a1 về đúng lưới t_ref (đề phòng solver adapt)
        a1u = interp1(t, a1, t_ref, 'linear', 'extrap');
        % RMSE theo SI (m/s^2)
        dif = a1u - a_ref;
        err = sqrt(mean(dif.^2));
        if any(~isfinite(err)) || err<=0, err = 1e6; end
        % Penalty nhẹ nếu overshoot quá dị (>80 g)
        if max(abs(a1u))/9.81 > 80, err = err + 1e4; end
    catch
        err = 1e7; % nếu ODE diverge
    end
end

function [t, a1] = sim_vehicle_accel(p, t_ref, v0, m1, m2)
% Chỉ trả về gia tốc xe để fitness nhanh
    [t, a1] = core_sim(p, t_ref, v0, m1, m2, true);
end

function [t, a1, a2, x1, x2, v1, v2] = sim_2dof_interp(p, t_ref, v0, m1, m2)
% Trả full outputs để vẽ
    [t, a1, a2, x1, x2, v1, v2] = core_sim(p, t_ref, v0, m1, m2, false);
end

function [tU, a1U, a2U, x1U, x2U, v1U, v2U] = core_sim(p, t_ref, v0, m1, m2, vehOnly)
% ODE45 solve rồi nội suy lên t_ref
    k1=p(1); c1=p(2); k2=p(3); c2=p(4);
    x0 = [0; v0; 0; v0];  % [x1 v1 x2 v2] tại t=0 cùng tốc độ v0
    odefun = @(tt,y) rhs(tt,y,m1,m2,k1,c1,k2,c2);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);
    [t, Y] = ode45(odefun, [t_ref(1) t_ref(end)], x0, opts);
    x1 = Y(:,1); v1 = Y(:,2);
    x2 = Y(:,3); v2 = Y(:,4);

    % Gia tốc từ phương trình động lực
    a1 = (-c1.*v1 - k1.*x1 + c2.*(v2 - v1) + k2.*(x2 - x1))./m1;
    a2 = (-c2.*(v2 - v1) - k2.*(x2 - x1))./m2;

    % Nội suy lên lưới chuẩn t_ref
    tU  = t_ref;
    a1U = interp1(t, a1, t_ref, 'linear','extrap');
    if vehOnly
        a2U=[]; x1U=[]; x2U=[]; v1U=[]; v2U=[];
        return;
    end
    a2U = interp1(t, a2, t_ref, 'linear','extrap');
    x1U = interp1(t, x1, t_ref, 'linear','extrap');
    x2U = interp1(t, x2, t_ref, 'linear','extrap');
    v1U = interp1(t, v1, t_ref, 'linear','extrap');
    v2U = interp1(t, v2, t_ref, 'linear','extrap');
end

function dy = rhs(~, y, m1, m2, k1, c1, k2, c2)
% 2-DOF: m1 với đất qua (k1,c1) nối tiếp, m2 nối với m1 qua (k2,c2)
    x1=y(1); v1=y(2); x2=y(3); v2=y(4);
    a1 = (-c1*v1 - k1*x1 + c2*(v2 - v1) + k2*(x2 - x1))/m1;
    a2 = (-c2*(v2 - v1) - k2*(x2 - x1))/m2;
    dy = [v1; a1; v2; a2];
end

function a = triangle_pulse(t, Apeak, Tp)
% Pulse mục tiêu dạng tam giác (đơn giản để demo)
    a = zeros(size(t));
    tt = t - t(1);
    up = tt<=Tp;
    a(up) = Apeak*(1 - tt(up)/Tp);
end
