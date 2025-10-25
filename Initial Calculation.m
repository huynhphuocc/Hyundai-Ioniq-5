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
%%
%II. Build the required calculation
%Specification
    %a. Dimension features
        L           = 4.635;        % [m]  Overall length
        W           = 1.890;        % [m]  Overall width
        H           = 1.600;        % [m]  Overall height
        GC          = 0.155;        % [m]  Ground clearance
        Wheelbase   = 2.999;        % [m]  Distance between front and rear axles
        CG_height   = 0.55;         % [m]  Estimated center of gravity height
        Weight      = 1905*9.81;    % [N]  Curb weight 
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
        Cabin_VolDensity = m_veh1 / Vol_Passenger;             % [kg/m^3] effective occupant space density
        Cargo_to_Passenger = Vol_Cargo_SeatsUp / Vol_Passenger;  % ratio indicator

    %e.Impact configuration:
        Impact_Velocity      = 64 / 3.6;                        % [m/s] initial impact speed (64 km/h)
        Impact_Angle         = 0;                               % [deg] head-on direction (0° deviation)
        Overlap_Ratio        = 0.40;                            % [-] 40% overlap between vehicle and barrier

    % Vehicle setup:
        v_veh1     = [Impact_Velocity, 0, 0];                   % [m/s] along X-axis
        Ground_Friction_Coeff= 0.80;                            % [-] coefficient of friction at tire–ground interface
        Suspension_Stiffness = 3.5e5;                           % [N/m] equivalent vertical stiffness (approx.)
        Suspension_Damping   = 2.0e3;                           % [N·s/m] equivalent damping

    % Simulation timing:
        End_time       = 0.2;                            % [s] estimated total impact duration
        Time_step      = 1e-5;                           % [s] analysis time increment (explicit solver reference)
