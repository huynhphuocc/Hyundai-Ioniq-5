clear ;
close all;
clc;
%{  
    PRELIMINARY CRASHWORTHINESS ESTIMATION - Hyundai Ioniq 5
    Purpose : Overall theoretical estimation before simulation
    Author  : Le Huynh Phuoc - 21145603 (Faculty of International Education – HCMUTE) 
%}
%% ------------------------------------------------------------------------------------------------

%   I.Initialize common format of figures/charts:
fig = figure('Name','PRELIMINARY CRASHWORTHINESS ESTIMATION - Hyundai Ioniq 5','NumberTitle','off');
    % Create a uitabgroup
    tabgroup = uitabgroup(fig);
    % Create uitab objects
    tab1 = uitab(tabgroup, 'Title', '1. Maxwell Model Data'); 
    tab2 = uitab(tabgroup, 'Title', '2. Deceleration');
    tab3 = uitab(tabgroup, 'Title', '3. Velocity');
    tab4 = uitab(tabgroup, 'Title', '4. Displacement');
    tab5 = uitab(tabgroup, 'Title', '5. Stiffness');
    tab6 = uitab(tabgroup, 'Title', '6. Damping');
    tab7 = uitab(tabgroup, 'Title', '7. Force');
    tab8 = uitab(tabgroup, 'Title', '8. Energy Summary');
    tab9 = uitab(tabgroup, 'Title', '9. Head Injury Metric');
    tab10= uitab(tabgroup, 'Title', '10. Chest Injury Metric');
    tab11= uitab(tabgroup, 'Title', '11. Neck Injury Metric');
    tab12= uitab(tabgroup, 'Title', '12. Femur Injury Metric');
    tab13= uitab(tabgroup, 'Title', '13. Vehicle Energy Distribution');
    tab14= uitab(tabgroup, 'Title', '14. Acceleration Severity Index');

%% -------------------------------
%II. Build the required calculation
%Specification
    %a. Dimension features
        Length      = 4.635;        % [m]  Overall length
        Width       = 1.890;        % [m]  Overall width
        Height      = 1.600;        % [m]  Overall height
        GC          = 0.155;        % [m]  Ground clearance
        Wheelbase   = 2.999;        % [m]  Distance between front and rear axles
        CG_height   = 0.55;         % [m]  Estimated center of gravity height
        g           = 9.81;         % [m/s] Gravitational acceleration
        Weight      = 1905*g;       % [N]  Curb weight 
        m_veh1      = 1905;         % [kg] Vehicle mass
        m_dum1      = 75;           % [kg] Dummy mass
        mi          = 0.70;         % Weight distribution coefficient (FWD)

    %b. Interior Dimensions (meter)
        Headroom_Front  = 1.011;    % [m] head room at front row
        Headroom_Rear   = 0.983;    % [m] head room at rear row
        Legroom_Front   = 1.059;    % [m] leg room at front row
        Legroom_Rear    = 1.001;    % [m] leg room at rear row
        Shoulder_Front  = 1.466;    % [m] shoulder room at front row
        Shoulder_Rear   = 1.466;    % [m] shoulder room at rear row
        Hip_Front       = 1.369;    % [m] hip room at front row
        Hip_Rear        = 1.362;    % [m] hip room at rear row

   %c. SAE Volume (cubic meter)
        Vol_Total           = 3.786;   % [m^3] total interior volume
        Vol_Passenger       = 3.016;   % [m^3] passenger compartment
        Vol_Cargo_SeatsUp   = 0.770;   % [m^3] rear cargo with seats up
        Vol_Cargo_SeatsDown = 1.679;   % [m^3] rear cargo with seats folded
        Vol_FrontStorage    = 0.024;   % [m^3] under-hood storage tray

   %d. Derived ratios (optional for future occupant/deformation study)
        Ratio_Headroom   = Headroom_Rear / Headroom_Front;   % rear/front comparison
        Ratio_Legroom    = Legroom_Rear / Legroom_Front;
        Cabin_VolDensity = m_veh1 / Vol_Passenger;               % [kg/m^3] effective occupant space density
        Cargo_to_Passenger = Vol_Cargo_SeatsUp / Vol_Passenger;  % ratio indicator

    %e.Impact configuration
        Impact_V1     = 50 /3.6;       % [m/s] frontal impact speed (50 km/h)
        Impact_V2     = 50 /3.6;       % [m/s] side impact speed (50 km/h)
        Impact_V3     = 32 /3.6;       % [m/s] side pole impact speed (50 km/h)
        Overlap       = 0.5;           % 50% overlap between vehicle and barrier

    %f. Vehicle setup
        Friction    = 0.80;       % coefficient of friction at tire–ground interface
        Crushmax    = 0.65;       % [m] Estimated maximum directional deformation
        eta         = 0.8;        % [%] Estimated efficiency of energy absorption
        v_rebound1  =-3.35;       % [m/s] Reference vehicle rebound speed 
     
    %g. Occupant setup
        v_rebound_OCC1   =-11.7;       % [m/s] Reference passenger rebound speed
        Occ_dmax         = 0.8;        % [m] Estimated torso displacement
        lamda            = 0.3;
        
    %h. Primary material inputs (HSS) 
        E   = 210e9;           % [Pa] Young's modulus (210 GPa)
        nu  = 0.30;            % Poisson (not used in 1D curve)
        rho = 7890;            % [kg/m^3] Density
        Sy  = 680e6;           % [Pa] Yield strength
        Sut = 760e6;           % [Pa] Ultimate tensile strength

    %i. Maxwell Model parameter
        K_coef  = (eta*m_veh1*Impact_V1^2)/Crushmax^2;      % [N/m] Stiffness coefficient
        C_coef  = lamda*0.5 *sqrt(K_coef*m_veh1);       % [Ns/m] Damping factor

    %k. Simulation timing
        End_time       = 0.2;                  % [s] estimated impact duration
        dt             = 0.001;                % [s] analysis time increment (1ms)
        t              = -0.001:dt:End_time;   % Range of time

%% -------------------------------
%1. APPROACH 1: Maxwell Model Data
%   Create axes, Switch to tab 1 using uitab(tab1);
    ax1 = axes('Parent', tab1);

%  Maxwell series' Initial condition 
x_maxwell = zeros(size(t)); v_maxwell = zeros(size(t)); a_maxwell = zeros(size(t));
v_maxwell(1) = Impact_V1;     

for i = 2:length(t) 
    a_maxwell(i)   = -(C_coef * v_maxwell(i-1) + K_coef*x_maxwell(i-1)) / m_veh1;       % F= Cv+Kx ; F= ma
    v_maxwell(i)   = v_maxwell(i-1) + a_maxwell(i)*dt;
    x_maxwell(i)   = x_maxwell(i-1) + v_maxwell(i)*dt;
end % Euler forward method

plot(t(2:end), a_maxwell(2:end)/g, 'LineWidth',1.5);
hold on
plot(t(2:end), v_maxwell(2:end)*3.6, '--');
plot(t(2:end), x_maxwell(2:end)*100);
ylabel('Deceleration (Gs) ; Velocity (km/h) ; Displacement (cm)','FontSize', 13);
xlabel('Duration (s)','FontSize', 13);
title('[Maxwell model] Crash Pulse - Velocity Reduction - Vehicle Displacement','FontSize', 16);
legend('Acceleration_{a}','Velocity_{v}','Displacement_{x}');
        % Find the highest point 
        [Max_a_lin, idx_Max_a_lin] = min(a_maxwell/9.81);
        Maximum_a_lin = [t(idx_Max_a_lin), Max_a_lin];
        % Find the end point
        Max_a_lin_end = [t(end), a_maxwell(end)/9.81];
        % Mark for 4 points
        h3=plot(Maximum_a_lin(1), Maximum_a_lin(2),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0]);
        set(get(get(h3,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        h4=plot(Max_a_lin_end(1), Max_a_lin_end(2),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
        set(get(get(h4,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        % Add text labels for the points
        text(Maximum_a_lin(1), Maximum_a_lin(2), sprintf('a_{peak} (%.1f Gs at %.f s)', Maximum_a_lin(2), Maximum_a_lin(1)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 11);
        % text(Max_a_lin_end(1), Max_a_lin_end(2), sprintf('a_{lin,end} (%.1f Gs at %.2f s)', Max_a_lin_end(2), Max_a_lin_end(1)), ...
        % 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 11);
        set(gca,'yMinorGrid','on');
        set(gca,'xMinorGrid','on');
        grid on;
        pbaspect([6 4 1]);
        hold off
%% 2. APPROACH 2: Genetic Algorithm Prediction Method  

% TARGET (from reference data)
t_target  = [0, 0.0894, End_time];
x1_target = [0, Crushmax, 0.2];
x2_target = [0, Occ_dmax, 0.4];
v1_target = [Impact_V1, 0, v_rebound1];
v2_target = [Impact_V1, 0, v_rebound_OCC1];

% Limitation of 26 parameters (ki1,ki2,ki3,ki4; xi1 xi2; ci1,ci2,ci3,ci4; vi1 vi2 Crushmax/OCC_dmax)
lb = [1e5,2e5,4e5,7.5e5,  0.15,0.45, 6e4,4e4,2e4,1e4,  5,10,  3e3,6e3,1e4,2e4, 0.02,0.3, 4e3,3e3,2e3,1e3, 1,8,   0.65,1 ]; % lower bound
ub = [2e5,4e5,7.5e5,1e6,  0.5,0.6,   1e5,6e4,4e4,2e4,  10,14, 6e3,1e4,2e4,5e4, 0.5,0.8,  2e4,4e3,3e3,2e3, 8,15,  0.65,1 ]; % upper bound

% Solving Genetic Algorithm of LPM
fprintf('Generating Genetic Algorithm ...\n');
% options = optimoptions('ga','PopulationSize',500,'MaxGenerations',1000,...
%     'EliteCount',30,'Display','iter','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-6,...
%     'UseParallel',true);
options = optimoptions('ga',...
    'PopulationSize',500,...
    'MaxGenerations',1500,...
    'MutationFcn',{@mutationadaptfeasible, 0.5},...   % seed mutation up to 0.5
    'CrossoverFraction',0.9,...                       
    'EliteCount',15,...
    'Display','iter',...
    'PlotFcn',@gaplotbestf,...
    'FunctionTolerance',1e-10,'UseParallel',true);

[best_params, best_err] = ga(@(p) cost_func(p,m_veh1,m_dum1,Impact_V1,t_target,x1_target,x2_target,v1_target,v2_target),...
    26,[],[],[],[],lb,ub,[],options);

% Solving ODE by Runge–Kutta methods
[t,y] = ode45(@(t,y) crash_ode(t,y,best_params,m_veh1,m_dum1,Impact_V1),[0 0.2],[0 Impact_V1 0 Impact_V1]);

% PLOT RESULT SIGNATURE
plot_results(t,y,best_params,Impact_V1,t_target,x1_target,x2_target,v1_target,v2_target, ...
    tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9,tab10,tab11,tab12,tab13,tab14);
fprintf('Completed... Error: %.2e\n', best_err);
fprintf('Best parameters:\n');
fprintf('%g ', best_params);
fprintf('\n');

% LOCAL FUNCTIONS 

function err = cost_func(p,m_veh1,m_dum1,Impact_V1,t_tgt,x1_tgt,x2_tgt,v1_tgt,v2_tgt)
    try
        [t,y] = ode45(@(t,y) crash_ode(t,y,p,m_veh1,m_dum1,Impact_V1),[0 0.2],[0 Impact_V1 0 Impact_V1]);
        if size(y,1)<50, err=1e12; return; end
    catch, err=1e12; return; end
    tq = 0:0.001:0.2;
    x1m = interp1(t,y(:,1),tq,'pchip',0); x2m = interp1(t,y(:,3),tq,'pchip',0);
    v1m = interp1(t,y(:,2),tq,'pchip',0); v2m = interp1(t,y(:,4),tq,'pchip',0);
    x1t = interp1(t_tgt,x1_tgt,tq,'pchip'); x2t = interp1(t_tgt,x2_tgt,tq,'pchip');
    v1t = interp1(t_tgt,v1_tgt,tq,'pchip'); v2t = interp1(t_tgt,v2_tgt,tq,'pchip');
    err = 5e8*mean((x1m-x1t).^2) + 5e8*mean((x2m-x2t).^2) + ...     % Mean squared error
          5e8*mean((v1m-v1t).^2) + 5e8*mean((v2m-v2t).^2);
    if max(x1m)>0.8 || max(x2m)>1.1, err = err + 1e10; end % Penalty if over-deformation
end

function dydt = crash_ode(~,y,p,m_veh1,m_dum1,~)
    x1=y(1); v1=y(2); x2=y(3); v2=y(4);
    k1 = piecewise_k1(x1,p(1:6),p(25));   % add Crushmax
    c1 = piecewise_c(v1,p(7:12));
    dx = x2-x1; dv = v2-v1;
    k2 = piecewise_k2(dx,p(13:18),p(26)); % add Occ_dmax
    c2 = piecewise_c(dv,p(19:24));
    Fstr  = k1*x1 + c1*v1;        % Deformation force
    Frest = k2*dx + c2*dv;        % Restraint system force
    dydt = [v1; (Frest-Fstr)/m_veh1; v2; -Frest/m_dum1]; % ODE 
end

function k = piecewise_k1(x,pk,Crushmax) % pk = [k1, k2, k3, k4, x1, x2]
    if x<=0, k=pk(1);
    elseif x<=pk(5), k=pk(1)+x*(pk(2)-pk(1))/pk(5);         
    elseif x<=pk(6), k=pk(2)+(x-pk(5))*(pk(3)-pk(2))/(pk(6)-pk(5));
    elseif x<=Crushmax, k=pk(3)+(x-pk(6))*(pk(4)-pk(3))/(Crushmax-pk(6));
    else, k=pk(4); end
end

function k = piecewise_k2(x,pk,Occ_dmax) % pk = [k1, k2, k3, k4, x1, x2]
    if x<=0, k=pk(1);
    elseif x<=pk(5), k=pk(1)+x*(pk(2)-pk(1))/pk(5);         
    elseif x<=pk(6), k=pk(2)+(x-pk(5))*(pk(3)-pk(2))/(pk(6)-pk(5));
    elseif x<=Occ_dmax, k=pk(3)+(x-pk(6))*(pk(4)-pk(3))/(Occ_dmax-pk(6));
    else, k=pk(4); end
end

function c = piecewise_c(v,pc)
    va = abs(v);
    v_th1= max(pc(5), 1e-3);              % v1
    v_th2= max(pc(6), v_th1 + 1e-3);      % v2 > v1
    v0= max(50/3.6, v_th2 + 1e-3);         % upper cap (≈impact speed)
    if va < 1e-3
        c=pc(1);
    elseif va <= v_th1
        c=pc(1)-(pc(1)-pc(2))/v_th1*va;
    elseif va <= v_th2
        c=pc(2)-(pc(2)-pc(3))/(v_th2-v_th1)*(va-v_th1);
    elseif va <= v0
        c=pc(3)-(pc(3)-pc(4))/(v0-v_th2)*(va-v_th2);   
    else
        c=pc(4);
    end
    c = max(c,1);
end

function plot_results(t,y,p,~,t_tgt,x1_tgt,x2_tgt,v1_tgt,v2_tgt, ...
    tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9,tab10,tab11,tab12,tab13,tab14)
    % Force all plots to render in main GUI figure
    fig = ancestor(tab2,'figure');
    figure(fig);
    tq = 0:0.001:0.2;
    x1m = interp1(t,y(:,1),tq); x2m = interp1(t,y(:,3),tq);
    v1m = interp1(t,y(:,2),tq); v2m = interp1(t,y(:,4),tq);
    % a1m = gradient(v1m, 0.001);
    % a2m = gradient(v2m, 0.001);
    a1m = gradient(smoothdata(v1m, 'gaussian', 20), 0.001);
    a2m = gradient(smoothdata(v2m, 'gaussian', 20), 0.001);
    kx1 = arrayfun(@(x) piecewise_k1(x,p(1:6),p(25))*x, x1m);
    cv1 = arrayfun(@(v) piecewise_c(v,p(7:12))*v, v1m);
    dx = x2m-x1m; dv = v2m-v1m;
    kx2 = arrayfun(@(d) piecewise_k2(d,p(13:18),p(26))*d, dx);
    cv2 = arrayfun(@(d) piecewise_c(d,p(19:24))*d, dv);
    if ~exist('results','dir'), mkdir('results'); end

    % Plot Acceleration curve
    ax2 = axes('Parent', tab2);    
    plot(ax2,tq,a1m/9.81,'LineWidth',1.5); 
    hold(ax2, 'on');
    plot(ax2,tq,a2m/9.81,'LineWidth',1.5); 
    ylabel(ax2,'Acceleration (Gs)','FontSize', 13);
    xlabel(ax2,'Duration (s)','FontSize', 13);
    legend(ax2,'a_1 (Vehicle)','a_2 (Passenger)'); 
    title(ax2,'[Genetic Algorithm Prediction] Crash Pulse','FontSize', 16);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax2, 'on');

    % Plot Velocity curve
    ax3 = axes('Parent', tab3);
    plot(ax3,tq,v1m,'LineWidth',1.5); hold(ax3,'on'); 
    plot(ax3,t_tgt,v1_tgt,'b--');  
    plot(ax3,tq,v2m,'LineWidth',1.5); 
    plot(ax3,t_tgt,v2_tgt,'r--'); 
    ylabel(ax3,'Velocity (m/s)','FontSize', 13);
    xlabel(ax3,'Duration (s)','FontSize', 13);
    legend(ax3,'v_1 (Vehicle)','v_1 (target)','v_2 (Passenger)','v_2 (target)'); 
    title(ax3,'[Genetic Algorithm Prediction] Velocity Reduction','FontSize', 16);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax3,'on');

    % Plot Displacement curve
    ax4 = axes('Parent', tab4);
    plot(ax4,tq,x1m,'LineWidth',1.5); 
    hold(ax4,'on'); 
    plot(ax4,t_tgt,x1_tgt,'b--'); 
    plot(ax4,tq,x2m,'LineWidth',1.5); 
    plot(ax4,t_tgt,x2_tgt,'r--'); 
    x_rel = x2m - x1m; 
    plot(ax4,tq,x_rel,'k-.','LineWidth',1.3);  % relative curve (torso)
    ylabel(ax4,'Displacement (m)','FontSize', 13);
    xlabel(ax4,'Duration (s)','FontSize', 13);
    legend(ax4,'x_1 (Vehicle)','x_1 (target)','x_2 (Passenger)','x_2 (target)', 'x_{torso}'); 
    title(ax4,'[Genetic Algorithm Prediction] Vehicle & Occupant Displacement','FontSize', 16);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax4,'on');

    % Plot Stiffness coefficient curve
    ax5 = axes('Parent', tab5);     
    plot(ax5,tq,kx1,'m','LineWidth',1.5); hold(ax5,'on'); 
    plot(ax5,tq,kx2,'c','LineWidth',1.5); 
    ylabel(ax5,'Stiffness Coefficient (N/m)','FontSize', 13);
    xlabel(ax5,'Duration (s)','FontSize', 13);
    legend(ax5,'K_{str}','K_{rest}'); 
    title(ax5,'[Genetic Algorithm Prediction] Stiffness Coefficient Characteristic','FontSize', 16);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax5,'on');

    % Plot Damping coefficient curve
    ax6 = axes('Parent', tab6);      
    plot(ax6,tq,abs(cv1),'m','LineWidth',1.5); hold(ax6,'on'); 
    plot(ax6,tq,abs(cv2),'c','LineWidth',1.5); 
    ylabel(ax6,'Damping Coefficient (Ns/m)','FontSize', 13);
    xlabel(ax6,'Duration (s)','FontSize', 13);
    legend(ax6,'C_{str}','C_{rest}'); 
    title(ax6,'[Genetic Algorithm Prediction] Damping Coefficient Characteristic','FontSize', 16);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax6,'on');

    % Plot Force curve
    ax7 = axes('Parent', tab7);      
    % F_str = kx1 + cv1;
    % F_rest = kx2 + cv2;
    F_str = arrayfun(@(x,v) piecewise_k1(x,p(1:6),p(25))*x + piecewise_c(v,p(7:12))*v, x1m, v1m);
    F_rest = arrayfun(@(dx,dv) piecewise_k2(dx,p(13:18),p(26))*dx + piecewise_c(dv,p(19:24))*dv, dx, dv);
    plot(ax7,tq, abs(F_str), 'LineWidth', 1.5); hold(ax7,'on');
    plot(ax7,tq, abs(F_rest), 'LineWidth', 1.5); 
    ylabel(ax7,'Force (N)','FontSize', 13);
    xlabel(ax7,'Duration (s)','FontSize', 13);
    legend(ax7,'Vehicle deformation force','Restraint system force'); 
    title(ax7,'[Genetic Algorithm Prediction] Force of Vehicle Deformation and Restraint System','FontSize', 14);
    set(gca,'yMinorGrid','on');
    set(gca,'xMinorGrid','on')
    grid(ax7,'on');  
end


%% ____ ENERGY SUMMARY_____________________ 
% VEHICLE kinetic energy
ax8 = axes('Parent', tab8);
KE_vehicle= 0.5*m_veh1 *(y(:,2).^2);
% OCCUPANT kinetic energy
KE_occ =0.5*m_dum1 *(y(:,4).^2);
% Recalculate stiffness and damping forces
x1 = y(:,1); v1 = y(:,2);
x2 = y(:,3); v2 = y(:,4);
dx = x2 - x1; dv = v2 - v1;
k1s = arrayfun(@(x) piecewise_k1(x, best_params(1:6), best_params(25)), x1);
c1s = arrayfun(@(v) piecewise_c( v, best_params(7:12)),v1);
k2s = arrayfun(@(d) piecewise_k2(d, best_params(13:18),best_params(26)), dx);
c2s = arrayfun(@(u) piecewise_c( u, best_params(19:24)),dv);
F_str  = k1s.*x1 + c1s.*v1;        % Vehicle deformation force
F_rest = k2s.*dx + c2s.*dv;        % Restraint system force
% Compute absorbed energy by integration of positive work (F*v)
P_str= F_str .* v1;              % Instantaneous power of structure
P_rest= F_rest .* dv;             % Instantaneous power of restraint
E_str= cumtrapz(t, max(P_str,0));   % Energy absorbed by vehicle structure
E_rest= cumtrapz(t, max(P_rest,0));   % Energy absorbed by restraint system
% Total energy balance
E_total = KE_vehicle + KE_occ + E_str + E_rest;
plot(ax8,t,KE_vehicle,'LineWidth',1.5); hold(ax8,'on');
plot(ax8,t,KE_occ,'LineWidth',1.5);
plot(ax8,t,E_str,'LineWidth',1.5);
plot(ax8,t,E_rest,'LineWidth',1.5);
ylabel(ax8,'Energy (J)','FontSize',13);
xlabel(ax8,'Duration (s)','FontSize',13);
legend(ax8,'Kinetic (Vehicle)','Kinetic (Occupant)',...
       'Absorbed (Vehicle)','Absorbed (Restraint)');
title(ax8,'[Energy Summary] Kinetic and Absorbed Energy','FontSize',16);
set(gca,'yMinorGrid','on');
set(gca,'xMinorGrid','on');
grid(ax8,'on');

%% 9. INJURY METRICS: HIC15 & HIC36 (based on a_2)
ax9 = axes('Parent', tab9);

% Acceleration of occupant
tq = t;                             % time vector from ODE
dt = mean(diff(tq));                % time step 
v2m = y(:,4);                       % occupant velocity (m/s)
a2m = gradient(v2m, dt);            % occupant accel (m/s^2)
a2g = abs(a2m/9.81);                % Gs (absolute)
%  Helper inline: calculate HIC with 'winSec' window 
IA = [0; cumsum(a2g(:))];           % prefix-sum (calculation following model)
N  = numel(a2g);
% --- HIC15 ---
win15   = 0.015;
Wmin15  = max(1, round(win15/dt));
HIC15   = 0; i1_15 = 1; i2_15 = Wmin15;

for i1 = 1:(N-Wmin15)
    i2max = min(N, i1 + ceil(1.5*Wmin15));   % to ~1.5*win
    for i2 = (i1+Wmin15):i2max
        dtw = tq(i2) - tq(i1);               % length of window (s)
        if dtw<=0 || dtw>win15, continue; end
        avg_a = (IA(i2)-IA(i1)) / (i2 - i1); % mean G
        hicv  = dtw * (avg_a^2.5);           % HIC formula
        if hicv > HIC15, HIC15 = hicv; i1_15 = i1; i2_15 = i2; end
    end
end
t1_15 = tq(i1_15); t2_15 = tq(i2_15);
% --- HIC36 ---
win36   = 0.036;
Wmin36  = max(1, round(win36/dt));
HIC36   = 0; i1_36 = 1; i2_36 = Wmin36;
for i1 = 1:(N-Wmin36)
    i2max = min(N, i1 + ceil(1.5*Wmin36));
    for i2 = (i1+Wmin36):i2max
        dtw = tq(i2) - tq(i1);
        if dtw<=0 || dtw>win36, continue; end
        avg_a = (IA(i2)-IA(i1)) / (i2 - i1); % mean (G)
        hicv  = dtw * (avg_a^2.5);
        if hicv > HIC36, HIC36 = hicv; i1_36 = i1; i2_36 = i2; end
    end
end
t1_36 = tq(i1_36); t2_36 = tq(i2_36);
% ===== Estimated AIS3+ & AIS level =====
approxRisk = @(hic) ( ...
    (hic<300).*0.05 + ...
    (hic>=300 & hic<700).*(0.05 + (hic-300)*(0.45/400)) + ...
    (hic>=700 & hic<1000).*(0.50 + (hic-700)*(0.40/300)) + ...
    (hic>=1000).*0.95 );
approxAIS  = @(hic) ( ...
    (hic<300).*1 + ...
    (hic>=300 & hic<500).*2 + ...
    (hic>=500 & hic<850).*3 + ...
    (hic>=850 & hic<1500).*4 + ...
    (hic>=1500).*5 );
risk15 = min(max(approxRisk(HIC15),0),1);
risk36 = min(max(approxRisk(HIC36),0),1);
ais15  = approxAIS(HIC15);
ais36  = approxAIS(HIC36);
    plot(ax9, tq, a2g, 'LineWidth', 1.5); hold(ax9,'on');
    ymax = max(a2g)*1.15 + eps;
% highlight HIC15 window
patch(ax9, [t1_15 t2_15 t2_15 t1_15], [0 0 ymax ymax], ...
      [0.95 0.75 0.75], 'FaceAlpha', 0.28, 'EdgeColor','none');
% highlight HIC36 window
patch(ax9, [t1_36 t2_36 t2_36 t1_36], [0 0 ymax ymax], ...
      [0.75 0.85 0.95], 'FaceAlpha', 0.28, 'EdgeColor','none');
% annotate edges 
plot(ax9, [t1_15 t2_15], [a2g(i1_15) a2g(i2_15)], 'r-', 'LineWidth', 1.0);
plot(ax9, [t1_36 t2_36], [a2g(i1_36) a2g(i2_36)], 'b-', 'LineWidth', 1.0);
txt15 = sprintf('HIC15 = %.1f | AIS~%d | Risk(AIS3+)~%d%%', HIC15, ais15, round(risk15*100));
txt36 = sprintf('HIC36 = %.1f | AIS~%d | Risk(AIS3+)~%d%%', HIC36, ais36, round(risk36*100));
text(t2_15, ymax*0.82, txt15, 'Color','r','FontSize',11,'HorizontalAlignment','right','Parent',ax9);
text(t2_36, ymax*0.68, txt36, 'Color','b','FontSize',11,'HorizontalAlignment','right','Parent',ax9);
    ylabel(ax9,'Occupant Head Acceleration (Gs)','FontSize', 13);
    xlabel(ax9,'Duration (s)','FontSize', 13);
    legend(ax9,'a_2 (Passenger)','HIC15 window','HIC36 window','Location','best');
    title(ax9,'[Injury Metric] HIC15 & HIC36','FontSize', 16);
    set(ax9,'yMinorGrid','on','xMinorGrid','on');
    grid(ax9,'on');
% Print to command window
fprintf('HIC15 = %.1f (%.3f–%.3f s), AIS~%d, Risk(AIS3+)~%d%%\n', HIC15, t1_15, t2_15, ais15, round(risk15*100));
fprintf('HIC36 = %.1f (%.3f–%.3f s), AIS~%d, Risk(AIS3+)~%d%%\n', HIC36, t1_36, t2_36, ais36, round(risk36*100));

%% 10. Chest Deflection & V*C
ax10 = axes('Parent', tab10);

x1m = interp1(t, y(:,1), tq, 'pchip');
x2m = interp1(t, y(:,3), tq, 'pchip');
% Chest deflection proxy = relative motion of torso wrt vehicle (positive compression)
x_rel = max(0, x2m - x1m);          % m
% Compression rate
v_rel = gradient(x_rel, dt);        % m/s
v_rel = max(0, v_rel);              % taking compression phase
% Thorax depth (50th male ~ 0.228 m). 
D_chest = 0.228;
try
    if evalin('base','exist(''Thorax_Depth'',''var'')'), D_chest = evalin('base','Thorax_Depth'); end
end
% Viscous Criterion
C   = x_rel / D_chest;              % compression ratio (−)
VC  = v_rel .* C;                   % m/s
[VCmax, idxVC] = max(VC);
tVC = tq(idxVC);
% Max chest deflection
[xmax, idxX] = max(x_rel);
tX = tq(idxX);
%  Plot with 2 y-axis 
yyaxis(ax10,'left');
    plot(ax10, tq, x_rel*1000, 'LineWidth', 1.5);           % mm
    yline(ax10, 55, '--', '55 mm limit', ...
          'LabelHorizontalAlignment','left','LabelVerticalAlignment','top');
    ylabel(ax10,'Chest deflection (mm)','FontSize',13);
    yyaxis(ax10,'right');
    plot(ax10, tq, VC, 'LineWidth', 1.5);                    % m/s
    yline(ax10, 1.0, '--', 'VC = 1.0 m/s', ...
          'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
    ylabel(ax10,'V*C (m/s)','FontSize',13);
% Markers + text
hold(ax10,'on');
    yyaxis(ax10,'left');  plot(ax10, tX, xmax*1000, 'ko','MarkerSize',5,'LineWidth',1);
    text(tX, xmax*1000, sprintf('  x_{max}=%.1f mm @ %.3fs', xmax*1000, tX), ...
         'FontSize',10,'Parent',ax10,'VerticalAlignment','bottom');
    
    yyaxis(ax10,'right'); plot(ax10, tVC, VCmax, 'ks','MarkerSize',5,'LineWidth',1);
    text(tVC, VCmax, sprintf('  VC_{max}=%.2f m/s @ %.3fs', VCmax, tVC), ...
         'FontSize',10,'Parent',ax10,'VerticalAlignment','bottom');
    xlabel(ax10,'Duration (s)','FontSize',13);
    title(ax10,'[Injury Metric] Chest Deflection & V*C','FontSize',16);
    legend(ax10, 'Deflection','VC','Location','best');
    set(ax10,'xMinorGrid','on','yMinorGrid','on'); grid(ax10,'on');
% console for report
fprintf('Chest deflection max = %.1f mm at %.3fs\n', xmax*1000, tX);
fprintf('VC max = %.2f m/s at %.3fs (Thorax depth = %.0f mm)\n', VCmax, tVC, D_chest*1000);

%% 11. Neck metrics (NIC & proxy Nij)
ax11 = axes('Parent', tab11);

v_torso = v2m;                       % torso velocity (m/s)
a_torso = a2m;                       % torso accel (m/s^2)
% --- Head surrogate by 1dof filter (head retard wrt to chest ~15 ms) ---
tau_neck = 0.015;                    % time constant ~ 15 ms
v_head = zeros(size(v_torso));
for i = 2:numel(v_head)
    v_head(i) = v_head(i-1) + dt*((v_torso(i) - v_head(i-1))/tau_neck);
end
a_head = gradient(v_head, dt);
% --- NIC (Boström): NIC = 0.2*a_rel + 0.5*v_rel^2 ---
% 0.2 m ~ distance head–T1 according to definition of NIC
v_rel = v_head - v_torso;
a_rel = a_head - a_torso;
NIC   = 0.2*a_rel + 0.5*(v_rel.^2);  % [m^2/s^2]
[NICmax, iNIC] = max(NIC);
tNIC = tq(iNIC);
% Raw Nij "proxy" from 2DOF (educational purpose) ---
m_head = 4.5;                        % head mass (kg) 50th male
Lneck  = 0.10;                       % length of neck - arm (m) approx
Fz     = m_head * a_head;            % axial neck force proxy (N)
Fshear = m_head * (a_head - a_torso);% shear due to relative sliding
My     = Fshear * Lneck;             % bending moment proxy (N·m)
% Reference Hybrid III 50th 
Fint_ten = 6806;     % N (tension)
Fint_comp= 6160;     % N (compression)
Mint_flex= 310;      % N·m
Mint_ext = 135;      % N·m
% 4 Nij combination (taking maximum)
N_te = max(Fz,0)/Fint_ten + max(abs(My),0)/Mint_ext;    % tension–extension
N_tf = max(Fz,0)/Fint_ten + max(abs(My),0)/Mint_flex;   % tension–flexion
N_ce = max(-Fz,0)/Fint_comp + max(abs(My),0)/Mint_ext;  % compression–extension
N_cf = max(-Fz,0)/Fint_comp + max(abs(My),0)/Mint_flex; % compression–flexion
Nij  = max([N_te; N_tf; N_ce; N_cf], [], 1);
[Nijmax, iNij] = max(Nij);
tNij = tq(iNij);
%  Plot NIC (left) & Nij (right) 
    yyaxis(ax11,'left');
    plot(ax11, tq, NIC, 'LineWidth', 1.5);
    yline(ax11, 15, '--', 'NIC=15 (guideline)', ...
        'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle');
    ylabel(ax11,'NIC (m^2/s^2)','FontSize',13);
    yyaxis(ax11,'right');
    plot(ax11, tq, Nij, 'LineWidth', 1.5);
    yline(ax11, 1.0, '--', 'Nij=1.0 limit', ...
        'LabelHorizontalAlignment','left');
    ylabel(ax11,'Nij (proxy)','FontSize',13);
% Mark & annotate peaks
yyaxis(ax11,'left');  hold(ax11,'on');
    plot(ax11, tNIC, NICmax, 'ko','MarkerSize',5,'LineWidth',1);
    text(tNIC, NICmax, sprintf('  NIC_{max}=%.1f @ %.3fs', NICmax, tNIC), ...
         'FontSize',10,'VerticalAlignment','bottom');
    yyaxis(ax11,'right');
    plot(ax11, tNij, Nijmax, 'ks','MarkerSize',5,'LineWidth',1);
    text(tNij, Nijmax, sprintf('  Nij_{max}=%.2f @ %.3fs', Nijmax, tNij), ...
         'FontSize',10,'VerticalAlignment','bottom');
    xlabel(ax11,'Duration (s)','FontSize',13);
    title(ax11,'[Neck Metrics] NIC & Nij (educational proxy)','FontSize',16);
    legend(ax11,'NIC','Nij','Location','best');
    set(ax11,'xMinorGrid','on','yMinorGrid','on'); grid(ax11,'on');
    % console for report
fprintf('NICmax = %.1f at %.3fs\n', NICmax, tNIC);
fprintf('Nij_max (proxy) = %.2f at %.3fs\n', Nijmax, tNij);

%% 12. Femur Load (proxy, each leg)
ax12 = axes('Parent', tab12);

x_rel = max(0, x2m - x1m);          % torso moves forward wrt vehicle (m)
v_rel = gradient(x_rel, dt);        % relative compressive velocity (m/s)
v_rel = max(0, v_rel);              % taking compression phase
% --- Param "knee bolster / fascia" (remember to take it to workspace so as to override)
% Knee_Clearance = 0.080;   % (m) -> daskboard (around 80 mm)
Knee_Clearance = Legroom_Front - 0.95; % 0.95 m is normal reactive distance
Knee_k         = 2.0e5;   % (N/m) ≈ 2 kN at 10 mm
Knee_c         = 5.0e3;   % (N·s/m) contact damping (for shock absorption)
try
    if evalin('base','exist(''Knee_Clearance'',''var'')'), Knee_Clearance = evalin('base','Knee_Clearance'); end
    if evalin('base','exist(''Knee_k'',''var'')'),Knee_k= evalin('base','Knee_k'); end
    if evalin('base','exist(''Knee_c'',''var'')'),Knee_c= evalin('base','Knee_c'); end
end
%  Contact model 1D: generate force when touching knee bolster
pen   = max(0, x_rel - Knee_Clearance);     % deflection (m)
F_knee_total = Knee_k.*pen + Knee_c.*v_rel; % total force from both thighs (N)
F_knee_total = max(0, F_knee_total);
% Divide equally on both sides (assuming symmetry)
F_femur_each = 0.5 * F_knee_total;
% Limiting display to avoid unreasonable spikes
F_femur_each = min(F_femur_each, 12e3);     % clamp 12 kN/leg
%  Reference threshold
thr_good   = 7.6e3;   % 7.6 kN (Euro NCAP high score)
thr_limit  = 10.0e3;  % 10 kN (FMVSS 208 high limit)
% --- Plot
plot(ax12, tq, F_femur_each/1000, 'LineWidth', 1.5); hold(ax12, 'on'); % kN
yline(ax12, thr_good/1000, '--', '7.6 kN (good)',  'LabelHorizontalAlignment','left');
yline(ax12, thr_limit/1000,'--', '10 kN (limit)',  'LabelHorizontalAlignment','left');
% Mark peak
[Fpk, ipk] = max(F_femur_each);
tpk = tq(ipk);
plot(ax12, tpk, Fpk/1000, 'ko', 'MarkerSize',5,'LineWidth',1);
text(tpk, Fpk/1000, sprintf('  Peak = %.1f kN @ %.3fs', Fpk/1000, tpk), ...
     'FontSize',10,'VerticalAlignment','bottom','Parent',ax12);
xlabel(ax12,'Duration (s)','FontSize',13);
ylabel(ax12,'Femur axial load (kN per leg)','FontSize',13);
title(ax12,'[Lower-Extremity] Femur Load (proxy via knee bolster contact)','FontSize',16);
legend(ax12,'Femur load (each leg)','Location','best');
set(ax12,'xMinorGrid','on','yMinorGrid','on'); grid(ax12,'on');
% console for report
fprintf('Femur peak (each leg) = %.1f kN at %.3fs (clearance=%.0f mm, k=%.1f kN/mm)\n', ...
    Fpk/1000, tpk, Knee_Clearance*1000, Knee_k/1e6);
%% 13. Vehicle Energy Distribution 
ax13 = axes('Parent', tab13);

v1m = y(:,2);                       % occupant velocity (m/s)
KE_vehicle = 0.5*m_veh1 .*(v1m.^2);
% Initial energy (E0)
if exist('Impact_V1','var') && ~isempty(Impact_V1)
    E0 = 0.5*(m_veh1+m_dum1)*(Impact_V1^2);
else
    v0 = max(v1m(1), v2m(1));
    E0 = 0.5*(m_veh1 + m_dum1)*(v0^2);
end
E0 = max(E0, eps);
% Calculated as proportion
E_abs_pct = 100*(E_str(end)/ E0);      % Absorption %
E_ke_pct = 100 *(KE_vehicle(end)/ E0); % KE remaining %
% Bar chart
vals = [E_abs_pct, E_ke_pct, max(0, 100 - (E_abs_pct + E_ke_pct))];
labels = {'Absorbed (Vehicle)', 'Remaining KE', 'Numerical loss'};
colors = [0.3 0.6 0.9; 0.9 0.6 0.3; 0.7 0.7 0.7];

barh(ax13, vals, 'FaceColor','flat'); 
for i = 1:length(vals)
    b = barh(ax13, i, vals(i), 'FaceColor', colors(i,:), 'BarWidth', 0.5);
    hold(ax13, 'on');
    text(vals(i)+2, i, sprintf('%.1f%%', vals(i)), 'FontSize',11, 'VerticalAlignment','middle');
end
set(ax13,'YTick',1:length(labels),'YTickLabel',labels,'FontSize',12);
xlabel(ax13,'Percentage of E_0 (%)','FontSize',13);
title(ax13,'[Energy Summary] Vehicle Energy Distribution at Final Time','FontSize',16);
xlim(ax13,[0 120]); grid(ax13,'on');
% Print cho report
fprintf('\n[Vehicle Energy Distribution]\n');
fprintf('  Absorbed (Vehicle): %.1f%%\n', E_abs_pct);
fprintf('  Remaining KE:       %.1f%%\n', E_ke_pct);
fprintf('  Numerical loss:     %.1f%%\n', max(0, 100 - (E_abs_pct + E_ke_pct)));

%% 14. Acceleration Severity Index (ASI)
ax14 = axes('Parent', tab14);

% --- Acceleration data (from occupant or CG)
a2m = gradient(y(:,4), dt); % acceleration of occupant (m/s^2)
% --- Sliding average filter (ISO 6487 delta = 0.05s)
delta = 0.05; 
N = round(delta / dt);
a_avg = movmean(a2m, N);
% --- ASI computation (longitudinal only)
g = 9.81;
a_ref = 12 * g;        % limit for frontal direction
ASI_t = abs(a_avg) / a_ref;
ASI_peak = max(ASI_t);
    plot(ax14, t, ASI_t, 'LineWidth', 1.6);
    hold(ax14,'on');
    yline(1.0,'g--','Class A (safe)');
    yline(1.4,'y--','Class B');
    yline(1.9,'r--','Class C (danger)');
    xlabel(ax14,'Time (s)','FontSize',13);
    ylabel(ax14,'ASI','FontSize',13);
    title(ax14,'[ECE R94 / ISO 6487] Acceleration Severity Index','FontSize',16);
    legend(ax14,'Instantaneous ASI','Location','best');
    set(ax14,'xMinorGrid','on','yMinorGrid','on');
    grid(ax14,'on');
% Console summary
fprintf('\n--- ASI ANALYSIS ---\n');
fprintf('Peak ASI = %.3f\n', ASI_peak);
if ASI_peak <= 1
    fprintf('→ Safety Class A (Low injury risk)\n');
elseif ASI_peak <= 1.4
    fprintf('→ Safety Class B (Moderate injury risk)\n');
elseif ASI_peak <= 1.9
    fprintf('→ Safety Class C (High injury risk)\n');
else
    fprintf('→ Above Class C — Severe impact!\n');
end

