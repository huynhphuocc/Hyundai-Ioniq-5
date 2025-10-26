clear ;
close all;
clc;
%{  
    INITIAL CRASHWORTHINESS ESTIMATION – Hyundai Ioniq 5
    Purpose : Quick theoretical estimation before simulation
    Author  : Le Huynh Phuoc - 21145603 (Faculty of International Education – HCMUTE) 
%}

%% ------------------------------------------------------------------------------------------------

%   I.Initialize common format of figures/charts:
fig = figure('Name','Initial Calculation of Hyundai Ioniq 5 Crashworthiness','NumberTitle','off');
    % Create a uitabgroup
    tabgroup = uitabgroup(fig);
    % Create uitab objects
    tab1 = uitab(tabgroup, 'Title', 'Estimated Crash Pulse'); 
    tab2 = uitab(tabgroup, 'Title', 'Estimated Velocity Reduction');
    tab3 = uitab(tabgroup, 'Title', 'Estimated Vehicle Displacement');
    tab4 = uitab(tabgroup, 'Title', 'Impact Force');
    tab5 = uitab(tabgroup, 'Title', 'Kinetic Energy');
    tab6 = uitab(tabgroup, 'Title', 'Energy Absorption');
    tab7 = uitab(tabgroup, 'Title', 'Energy Summary');
    tab8 = uitab(tabgroup, 'Title', 'Estimated Deformation');
    tab9 = uitab(tabgroup, 'Title', 'Stress - Strain Curve'); 
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
        Impact_V2     = 60 /3.6;       % [m/s] side impact speed (50 km/h)
        Impact_V3     = 32 /3.6;       % [m/s] side pole impact speed (50 km/h)
        Overlap       = 0.5;           % 50% overlap between vehicle and barrier

    % Vehicle setup
        Ground_Friction_Coeff= 0.80;           % [-] coefficient of friction at tire–ground interface
        Suspension_Stiffness = 3.5e5;          % [N/m] equivalent vertical stiffness (approx.)
        Suspension_Damping   = 2.0e3;          % [N·s/m] equivalent damping

    % Maxwell Model parameter
        beta = 1.8;
        Ap_desired = 30*g;
        omega_n = Ap_desired/(beta*Impact_V1);
        K_coef  = omega_n^2 * m_veh1;            % [-] coefficient of friction at tire–ground interface
        C_coef  = 2 * omega_n * m_veh1;          % [N/m] equivalent vertical stiffness (approx.)
        K_eff = omega_n^2 * m_veh1;              % For Nonlinear model setup
        C_eff = 2 * omega_n * m_veh1;            % For Nonlinear model setup
    % Simulation timing
        End_time       = 0.4;                  % [s] estimated impact duration
        dt             = 0.001;                % [s] analysis time increment (1ms)
        t              = -0.001:dt:End_time;        % Range of time

%% -------------------------------
%1. Estimated Crash Pulse
%   Create axes, Switch to tab 1 using uitab(tab1);
    ax1 = axes('Parent', tab1);

% --- Maxwell series' Initial condition ---
x_lin = zeros(size(t)); v_lin = zeros(size(t)); a_lin = zeros(size(t));
x_non = zeros(size(t)); v_non = zeros(size(t)); a_non = zeros(size(t));
v_lin(1) = Impact_V1;     v_non(1) = Impact_V1;

Nk = 0.35;     % Nonlinear coefficient of k(x)
Nc = 0.01;     % Nonlinear coefficient of c(v)

for i = 2:length(t)
    % ===== Linear Maxwell series =====
    a_lin(i)   = -(C_coef * v_lin(i-1) + K_coef*x_lin(i-1)) / m_veh1;       % F= Cv+Kx ; F= ma
    v_lin(i)   = v_lin(i-1) + a_lin(i)*dt;
    x_lin(i)   = x_lin(i-1) + v_lin(i)*dt;
end
for i = 2:length(t)
    % ===== Nonlinear Maxwell series =====
    K_non = K_eff * (1 + Nk * x_non(i-1)^2);
    C_non = C_eff * (1 + Nc * abs(v_non(i-1)));
    
    a_non(i) = -(C_non/m_veh1)*v_non(i-1) - (K_non/m_veh1)*x_non(i-1);
    v_non(i) = v_non(i-1) + a_non(i)*dt;
    x_non(i) = x_non(i-1) + v_non(i)*dt;
end

plot(ax1, t(2:end), a_lin(2:end)/9.81, 'LineWidth',1.6);
hold on
plot(ax1, t(2:end), a_non(2:end)/9.81, '--', 'LineWidth',1.6, 'Color',[0.85 0.33 0.10]);
ylabel(ax1, 'Deceleration (Gs)','FontSize', 13);
xlabel(ax1, 'Duration (s)','FontSize', 13);
title(ax1, 'Estimated Crash Pulse - Linear vs Nonlinear Maxwell Model','FontSize', 14);
legend(ax1,'a_{linear}','v_{linear}','v_{nonlinear}','a_{nonlinear}');
        % Find the highest point 
        [Max_a_lin, idx_Max_a_lin] = min(a_lin/9.81);
        Maximum_a_lin = [t(idx_Max_a_lin), Max_a_lin];
        [Max_a_non, idx_Max_a_non] = min(a_non/9.81);
        Maximum_a_non = [t(idx_Max_a_non), Max_a_non];
        % Find the end point
        Max_a_lin_end = [t(end), a_lin(end)/9.81];
        Max_a_non_end = [t(end), a_non(end)/9.81];
        % Mark for 4 points
        h3=plot(Maximum_a_lin(1), Maximum_a_lin(2),'o','MarkerSize',5,'MarkerFaceColor',[1 0 0]);
        set(get(get(h3,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        h4=plot(Max_a_lin_end(1), Max_a_lin_end(2),'o','MarkerSize',5,'MarkerFaceColor',[0 0 0]);
        set(get(get(h4,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        h5=plot(Maximum_a_non(1), Maximum_a_non(2),'^','MarkerSize',5,'MarkerFaceColor',[0.85 0.33 0.10]);
        set(get(get(h5,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        h6=plot(Max_a_non_end(1), Max_a_non_end(2),'^','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5]);
        set(get(get(h6,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        % Add text labels for the points
        text(Maximum_a_lin(1), Maximum_a_lin(2), sprintf('a_{lin,peak} (%.1f Gs at %.f s)', Maximum_a_lin(2), Maximum_a_lin(1)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 11);
        text(Max_a_lin_end(1), Max_a_lin_end(2), sprintf('a_{lin,end} (%.1f Gs at %.2f s)', Max_a_lin_end(2), Max_a_lin_end(1)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 11);
        text(Maximum_a_non(1), Maximum_a_non(2), sprintf('a_{non,peak} (%.1f Gs at %.f s)', Maximum_a_non(2), Maximum_a_non(1)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 11);
        text(Max_a_non_end(1), Max_a_non_end(2), sprintf('a_{non,end} (%.1f Gs at %.2f s)', Max_a_non_end(2), Max_a_non_end(1)), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 11);
        set(gca,'yMinorGrid','on');
        set(gca,'xMinorGrid','on');
        grid on;
        pbaspect([6 4 1]);
        hold(ax1, 'off');

%% -------------------------------
%2. Estimated Velocity Reduction
    ax2 = axes('Parent', tab2);

plot(t(2:end), v_lin(2:end), 'LineWidth',1.6);
hold on
plot(t(2:end), v_non(2:end), '--', 'LineWidth',1.6);
ylabel(ax1, 'Deceleration (Gs)');
xlabel(ax1, 'Duration (s)');
title(ax1, 'Estimated Crash Pulse - Linear vs Nonlinear Maxwell Model');
legend(ax1,'a_{linear}','v_{linear}','v_{nonlinear}','a_{nonlinear}');
        set(gca,'yMinorGrid','on');
        set(gca,'xMinorGrid','on');
        grid on;
        pbaspect([6 4 1]);
        hold off;
