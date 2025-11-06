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
    tab9 = uitab(tabgroup, 'Title', '9. Stress - Strain Curve'); 
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
        Occ_dmax         = 1;          % [m] Estimated passenger displacement
        lamda            = 0.3;
        
    %h. Primary material inputs (HSS) 
        E   = 210e9;           % [Pa] Young's modulus (210 GPa)
        nu  = 0.30;            % Poisson (not used in 1D curve)
        rho = 7890;            % [kg/m^3] Density
        Sy  = 680e6;           % [Pa] Yield strength
        Sut = 760e6;           % [Pa] Ultimate tensile strength

    %i. Maxwell Model parameter
        K_coef  = (eta*m_veh1*Impact_V1^2)/Crushmax^2;      % [N/m] Stiffness coefficient
        C_coef  = lamda* 0.5 * sqrt(K_coef * m_veh1);       % [Ns/m] Damping factor

    %k. Simulation timing
        End_time       = 0.2;                  % [s] estimated impact duration
        dt             = 0.001;                % [s] analysis time increment (1ms)
        t              = -0.001:dt:End_time;        % Range of time

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
x1_target = [0, Crushmax, 0];
x2_target = [0, Occ_dmax, 0];
v1_target = [Impact_V1, 0, v_rebound1];
v2_target = [Impact_V1, 0, v_rebound_OCC1];

% Limitation of 24 parameters (ki1,ki2,ki3,ki4; ci1,ci2,ci3,ci4; xi1 xi2;
% vi1 vi2)
lb = [1e6,5e5,1e5,1e4,0.1,0.3,0.7269,5e4,1e4,5e3,500,2,1e6,5e5,1e5,1e4,0.03,0.1,0.2731,5e3,1e3,500,50,1];
ub = [1e8,5e7,1e7,1e6,0.4,0.6,0.7269,5e6,1e6,5e5,5e4,10,1e8,5e7,1e7,1e6,0.15,0.2,0.2731,5e5,1e5,5e4,5e3,8];

% Solving Genetic Algorithm of LPM
fprintf('Generating Genetic Algorithm ...\n');
options = optimoptions('ga','PopulationSize',300,'MaxGenerations',1000,...
    'EliteCount',30,'Display','iter','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-6);
[best_params, best_err] = ga(@(p) cost_func(p,m_veh1,m_dum1,Impact_V1,t_target,x1_target,x2_target,v1_target,v2_target),...
    24,[],[],[],[],lb,ub,[],options);

% Solving ODE by Runge–Kutta methods
[t,y] = ode45(@(t,y) crash_ode(t,y,best_params,m_veh1,m_dum1,Impact_V1),[0 0.2],[0 Impact_V1 0 Impact_V1]);

% PLOT RESULT
plot_results(t,y,best_params,Impact_V1,t_target,x1_target,x2_target,v1_target,v2_target,tab2,tab3,tab4,tab5,tab6,tab7);
fprintf('Completed... Error: %.2e\n', best_err);


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
    err = 1e8*mean((x1m-x1t).^2) + 1e8*mean((x2m-x2t).^2) + ...
          1e9*mean((v1m-v1t).^2) + 1e9*mean((v2m-v2t).^2);
    if max(x1m)>0.8 || max(x2m)>1.1, err = err + 1e10; end
end

function dydt = crash_ode(~,y,p,m_veh1,m_dum1,~)
    x1=y(1); v1=y(2); x2=y(3); v2=y(4);
    k1 = piecewise_k(x1,p(1:7));
    c1 = piecewise_c(v1,p(8:12));
    dx = x2-x1; dv = v2-v1;
    k2 = piecewise_k(dx,p(13:19));
    c2 = piecewise_c(dv,p(20:24));
    Fstr  = k1*x1 + c1*v1;        
    Frest = k2*dx + c2*dv;        
    dydt = [v1; (Frest-Fstr)/m_veh1; v2; -Frest/m_dum1];  
end

function k = piecewise_k(x,pk)
    if x<=0, k=pk(1);
    elseif x<=pk(5), k=pk(1)+(pk(2)-pk(1))/pk(5)*x;
    elseif x<=pk(6), k=pk(2)+(pk(3)-pk(2))/(pk(6)-pk(5))*(x-pk(5));
    elseif x<=pk(7), k=pk(3)+(pk(4)-pk(3))/(pk(7)-pk(6))*(x-pk(6));
    else, k=pk(4); end
end

function c = piecewise_c(v,pc)
    va = abs(v);
    if va<1e-3, c=pc(1);
    elseif va<=pc(5), c=pc(1)-(pc(1)-pc(2))/pc(5)*va;
    elseif va<=10, c=pc(2)-(pc(2)-pc(3))/(10-pc(5))*(va-pc(5));
    else, c=pc(3)-(pc(3)-pc(4))/(13.89-10)*(va-10); end
    c = max(c,1);
end

function plot_results(t,y,p,~,t_tgt,x1_tgt,x2_tgt,v1_tgt,v2_tgt,tab2,tab3,tab4,tab5,tab6,tab7)
    % Force all plots to render in main GUI figure
    fig = ancestor(tab2,'figure');
    figure(fig);

    tq = 0:0.001:0.2;
    x1m = interp1(t,y(:,1),tq); x2m = interp1(t,y(:,3),tq);
    v1m = interp1(t,y(:,2),tq); v2m = interp1(t,y(:,4),tq);
    a1m = gradient(v1m, 0.001);
    a2m = gradient(v2m, 0.001);
    kx1 = arrayfun(@(x) piecewise_k(x,p(1:7))*x, x1m);
    cv1 = arrayfun(@(v) piecewise_c(v,p(8:12))*v, v1m);
    dx = x2m-x1m; dv = v2m-v1m;
    kx2 = arrayfun(@(d) piecewise_k(d,p(13:19))*d, dx);
    cv2 = arrayfun(@(d) piecewise_c(d,p(20:24))*d, dv);
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
    grid(ax2, 'on');

    % Plot Velocity curve
    ax3 = axes('Parent', tab3);
    plot(ax3,tq,v1m,'LineWidth',1.5); hold(ax3,'on'); 
    plot(ax3,t_tgt,v1_tgt,'b--');  
    plot(ax3,tq,v2m,'LineWidth',1.5); 
    plot(ax3,t_tgt,v2_tgt,'b--'); 
    ylabel(ax3,'Velocity (m/s)','FontSize', 13);
    xlabel(ax3,'Duration (s)','FontSize', 13);
    legend(ax3,'v_1 (Vehicle)','v_1 (target)','v_2 (Passenger)','v_2 (target)'); 
    title(ax3,'[Genetic Algorithm Prediction] Velocity Reduction','FontSize', 16);
    grid(ax3,'on');

    % Plot Displacement curve
    ax4 = axes('Parent', tab4);
    plot(ax4,tq,x1m,'LineWidth',1.5); 
    hold(ax4,'on'); 
    plot(ax4,t_tgt,x1_tgt,'b--'); 
    plot(ax4,tq,x2m,'LineWidth',1.5); 
    plot(ax4,t_tgt,x2_tgt,'b--'); 
    ylabel(ax4,'Displacement (m)','FontSize', 13);
    xlabel(ax4,'Duration (s)','FontSize', 13);
    legend(ax4,'x_1 (Vehicle)','x_1 (target)','x_2 (Passenger)','x_2 (target)'); 
    title(ax4,'[Genetic Algorithm Prediction] Vehicle & Occupant Displacement','FontSize', 16);
    grid(ax4,'on');

    % Plot Stiffness coefficient curve
    ax5 = axes('Parent', tab5);     
    plot(ax5,tq,kx1,'m','LineWidth',1.5); hold(ax5,'on'); 
    plot(ax5,tq,kx2,'c','LineWidth',1.5); 
    ylabel(ax5,'Stiffness Coefficient (N/m)','FontSize', 13);
    xlabel(ax5,'Duration (s)','FontSize', 13);
    legend(ax5,'F_{str,k}','F_{rest,k}'); 
    title(ax5,'[Genetic Algorithm Prediction] Stiffness Coefficient Characteristic','FontSize', 16);
    grid(ax5,'on');

    % Plot Damping coefficient curve
    ax6 = axes('Parent', tab6);      
    plot(ax6,tq,cv1,'m','LineWidth',1.5); hold(ax6,'on'); 
    plot(ax6,tq,cv2,'c','LineWidth',1.5); 
    ylabel(ax6,'Damping Coefficient (Ns/m)','FontSize', 13);
    xlabel(ax6,'Duration (s)','FontSize', 13);
    legend(ax6,'F_{str,c}','F_{rest,c}'); 
    title(ax6,'[Genetic Algorithm Prediction] Damping Coefficient Characteristic','FontSize', 16);
    grid(ax6,'on');

    % Plot Force curve
    ax7 = axes('Parent', tab7);      
    F_str = kx1 + cv1;
    F_rest = kx2 + cv2;
    plot(ax7,tq, F_str, 'LineWidth', 1.5); hold(ax7,'on');
    plot(ax7,tq, F_rest, 'LineWidth', 1.5); 
    ylabel(ax7,'Force (N)','FontSize', 13);
    xlabel(ax7,'Duration (s)','FontSize', 13);
    legend(ax7,'Vehicle deformation force','Restraint system force'); 
    title(ax7,'[Genetic Algorithm Prediction] Force of Vehicle Deformation and Restraint System','FontSize', 14);
    grid(ax7,'on');
end

