clc; clear all; close all;
%% PART TWO: Time Optimal and Range Optimal Aircraft Launches of SpaceShipTwo and GOLauncher One

% The time optimal simulation portion of this script is derived from Appendix B from "Optimal Controls With
% Aerospace Applications". 
% Drag on the vehicle and varying mass are calculated for time optimal trajectory. 

%% Time-Optimal Trajectory to LEO for SpaceShipTwo and GOLauncher
for iter = 1:2
    global g Vc h K1 scale F m0 mdot xbar0 ybar0 Vxbar0 Vybar0 ybarf Vxbarf Vybarf
    % Boundary Conditions for SpaceShipTwo and GOLauncher 
    %SpaceShipTwo Conditions
    if iter == 1 
        h_init = 15000; % m ,initial altitude
        V_init = 55; % m/s, initial velocity
        h = 150000; % m, final altitude 
        Vc = sqrt(3.9860044e5/(6378.14+h/1000))*1000; % m/s, final circular speed
        g = 9.80665; % m/s^2, gravity
        F = 2.7e6; % N, constant thrust
        h_scale = 8440; % m, atmospheric scale-height 
        scale = h/h_scale; % constant to simplify EOM
        rho_ref = 1.225; % reference density
        r = 2.3/2; %m, radius of crew cabin
        A = pi*r^2; % m^22, cross-sectional area

        % Initial conditions
        xbar0 = 0; % initial x-position
        ybar0 = h_init/h; % initial y-position
        Vxbar0 = V_init/Vc; % initial downrange velocity
        Vybar0 = 0; % initial vertical velocity

        % Final conditions
        ybarf = h/h; % final altitude
        Vxbarf = Vc/Vc; % final downrange velocity
        Vybarf = 0; % final vertical velocity
        Hf = -1; %final value of Hamiltonian
    end
    %GOLauncher Conditions
    if iter == 2
        h_init = 10e3; % m ,initial altitude
        V_init = 220; % m/s, initial velocity
        h = 150000; % m, final altitude 
        Vc = sqrt(3.9860044e5/(6378.14+h/1000))*1000; % m/s, final circular speed
        g = 9.80665; % m/s^2, gravity
        F = 2.5e6; % N, constant thrust
        h_scale = 8440; % m, atmospheric scale-height 
        scale = h/h_scale; % constant to simplify EOM
        rho_ref = 1.225; % reference density
        r = .5/2; %m, radius of crew cabin
        A = pi*r^2; % m^22, cross-sectional area

        % Initial conditions
        xbar0 = 0; % initial x-position
        ybar0 = h_init/h; % initial y-position
        Vxbar0 = V_init/Vc; % initial downrange velocity
        Vybar0 = 0; % initial vertical velocity

        % Final conditions
        ybarf = h/h; % final altitude
        Vxbarf = Vc/Vc; % final downrange velocity
        Vybarf = 0; % final vertical velocity
        Hf = -1; %final value of Hamiltonian
    end
  

% Solution
    % CONSTANT MASS / NO DRAG Parameters
    m0 = 60880; % kg, average mass of rocket
    CD = 0; % no drag case has Cd = 0
    mdot = 0;
    K1 = rho_ref*CD*A*Vc/(2*m0); %constant to simplify EOM 

    % Initial Guesses
    t0 = 0;
    yinit = [xbar0 ybar0 Vxbar0 Vybar0 0 -1 0]; %lambda20 lambda30 lambda40
    tf_guess = 500; % sec, initial guess for final time
    Nt = 80;
    tau = linspace(0,1,Nt)'; % nondimensional time vector

    solinit = bvpinit(tau,yinit,tf_guess);

    % Solution - CONSTANT MASS / NO DRAG

    sol = bvp4c(@ascent_odes_tf, @ascent_bcs_tf, solinit);
    % Extract the final time from the solution:
    tf = sol.parameters(1);
    % Evaluate the solution 
    Z = deval(sol,tau);
    % Convert back to dimensional time for plotting
    time = t0 + tau.*(tf-t0);

    x_sol = Z(1,:)*h/1000;
    y_sol = Z(2,:)*h/1000;
    vx_sol = Z(3,:)*Vc/1000;
    vy_sol = Z(4,:)*Vc/1000;
    lambda2_bar_sol = Z(5,:);
    lambda3_bar_sol = Z(6,:);
    lambda4_bar_sol = Z(7,:);

    % Solution - VARIABLE MASS / NO DRAG
    if iter == 1
        m0 = 13154.179;
        mdot = 112.65;
        delta_tf = 114;
        CD = 0; 
        K1 = rho_ref*CD*A*Vc/(2*m0);
    end
    if iter == 2
        m0 = 5700;
        mdot = 65;
        delta_tf = 154;
        CD = 0; 
        K1 = rho_ref*CD*A*Vc/(2*m0);
    end

    solinit_mass = solinit;
    solinit_mass.y = Z;
    tf = sol.parameters(1);
    solinit_mass.parameters(1) = tf-delta_tf;

    sol_mass = bvp4c(@ascent_odes_tf, @ascent_bcs_tf, solinit_mass);
    % Extract the final time from the solution:
    tf_mass = sol_mass.parameters(1);
    % Evaluate the solution 
    Z_mass = deval(sol_mass,tau);
    % Convert back to dimensional time for plotting
    time_mass = t0 + tau.*(tf-t0);

    x_sol_mass = Z_mass(1,:)*h/1000;
    y_sol_mass = Z_mass(2,:)*h/1000;
    vx_sol_mass = Z_mass(3,:)*Vc/1000;
    vy_sol_mass = Z_mass(4,:)*Vc/1000;
    lambda2_bar_sol_mass = Z_mass(5,:);
    lambda3_bar_sol_mass = Z_mass(6,:);
    lambda4_bar_sol_mass = Z_mass(7,:);


    % Solution - VARIABLE MASS / WITH DRAG
    if iter == 1
        m0 = 13154.179;
        mdot = 112.65;
        delta_tf = 114;
        CD = .3; 
        K1 = rho_ref*CD*A*Vc/(2*m0);
    end
    if iter == 2
        m0 = 5700;
        mdot = 65;
        delta_tf = 154;
        CD = .1; 
        K1 = rho_ref*CD*A*Vc/(2*m0);
    end

    solinit_mass_drag = solinit_mass;
    solinit_mass_drag.y = Z_mass;
    tf_mass = sol_mass.parameters(1);
    solinit_mass_drag.parameters(1) = tf_mass;

    sol_mass_drag = bvp4c(@ascent_odes_tf, @ascent_bcs_tf, solinit_mass_drag);
    % Extract the final time from the solution:
    tf_mass_drag = sol_mass_drag.parameters(1);
    % Evaluate the solution 
    Z_mass_drag = deval(sol_mass_drag,tau);
    % Convert back to dimensional time for plotting
    time_mass_drag = t0 + tau.*(tf-t0);

    x_sol_mass_drag = Z_mass_drag(1,:)*h/1000;
    y_sol_mass_drag = Z_mass_drag(2,:)*h/1000;
    vx_sol_mass_drag = Z_mass_drag(3,:)*Vc/1000;
    vy_sol_mass_drag = Z_mass_drag(4,:)*Vc/1000;
    lambda2_bar_sol_mass_drag = Z_mass_drag(5,:);
    lambda3_bar_sol_mass_drag = Z_mass_drag(6,:);
    lambda4_bar_sol_mass_drag = Z_mass_drag(7,:);

% Plots 
    figure(1)
    subplot(2,2,1); hold on
    plot(time_mass_drag,x_sol_mass_drag); grid on
    xlabel('Time (sec)')
    ylabel('x (km)')
    %xlim([t0 tf])
    legend('SpaceShipTwo', 'GOLauncher', 'location', 'northwest')

    subplot(2,2,2); hold on
    plot(time_mass_drag,y_sol_mass_drag); grid on
    xlabel('Time (sec)')
    ylabel('y (km)')
    %xlim([t0 tf])
    legend('SpaceShipTwo', 'GOLauncher', 'location', 'northwest')

    subplot(2,2,3); hold on
    plot(time_mass_drag,vx_sol_mass_drag); grid on
    xlabel('Time (sec)')
    ylabel('V_x (km/s)')
    %xlim([t0 tf])
    legend('SpaceShipTwo', 'GOLauncher', 'location', 'northwest')

    subplot(2,2,4); hold on
    plot(time_mass_drag,vy_sol_mass_drag); grid on
    xlabel('Time (sec)')
    ylabel('V_y (km/s)')
    %xlim([t0 tf])
    legend('SpaceShipTwo', 'GOLauncher', 'location', 'northwest')

    figure(2); hold on; grid on
    plot(time_mass_drag,atand(lambda4_bar_sol_mass_drag./lambda3_bar_sol_mass_drag)); grid on
    xlabel('Time (sec)')
    ylabel('\theta (deg)')
    title('Steering Angle - Dante Sanaei');
    %xlim([t0 tf])
    legend('SpaceShipTwo', 'GOLauncher')

    figure(3); hold on; grid on
    plot(x_sol_mass_drag,y_sol_mass_drag);
    xlabel('Downrange Position, x (km)');
    ylabel('Altitude, y (km)');
    title('Altitude vs. Downrange Position - Dante Sanaei');
    %xlim([x_sol_mass_drag(1) x_sol_mass_drag(end)]);
    ylim([0 200]);
    legend('SpaceShipTwo', 'GOLauncher')

% Final Parameters
    final_x = x_sol_mass_drag(end)
    final_y = y_sol_mass_drag(end)
    final_vx = vx_sol_mass_drag(end)
    final_vy = vy_sol_mass_drag(end)
    final_time = tf
    
    
    clear all
end

%% Optimal Range Trajectory For SpaceShipTwo and GOLauncher
time_in_space_S1 = [];
time_in_space_GO = [];
angles_S1 = [];
angles_GO = [];

for iter = 1:2
    if iter == 1
        g = 9.81;
        mdot = 112;
        F = 270e3;
        Ispm = F/mdot;
        m = 13154.179/1.35;                                                           
        f = (mdot * Ispm)/m
        T = 87;
        g2f = g/f;
        y_init = 15e3;
        V0 = 55;
    end
    if iter == 2
        g = 9.81;
        mdot = 65;
        F = 250e3;
        Ispm = F/mdot;
        m = 5700/1.35;
        f = (mdot * Ispm)/m;
        T = 60;
        g2f = g/f;
        y_init = 10e3;
        V0 = 220;
    end
    
    for i = 0:.001:pi/2
        anglform = g2f*sin(i)^3 - 2*sin(i)^2 + 1;
        i;
        if  anglform <= 0.001 && anglform >= -.001
            anglform;
            optimal_range_theta = i;
            rad2deg(i);
        end
    end        
    final_time = 800;

    for angle = 0:1:80
        theta = deg2rad(angle);
        Vx0 = V0*cos(theta);
        Vy0 = V0*sin(theta);
        Vx1 = (f*T*cos(theta)) + Vx0;
        Vy1 = (f*sin(theta)-g)*T + Vy0;
        x1 = .5*f*T^2*cos(theta);
        y1 = .5*(f*sin(theta)-g)*T^2 + y_init;

        if  isreal(-(sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) - Vy1) / g)
            u = sqrt(Vx1^2+Vy1^2);
            x_burn = 0:.1:x1;
            slope = (y1-y_init) / x1;
            y_burn = slope * x_burn + y_init;
            time = 0:.1:final_time;
            x_coast = x1 + Vx1*time ;
            y_coast = y1 + Vy1.*time - .5 * g * time.^2;
            k = find(y_coast >-.01,1, 'last');
       
            x = [x_burn x_coast(1:k)] / 1000;
            y = [y_burn y_coast(1:k)] / 1000;
            
            if iter == 1
                figure(4)
                plot(x,y); grid on; hold on;
                xlabel('Downrange Position, x (km)');
                ylabel('Altitude, y (km)');
                title('SpaceShipTwo Ballistic Trajectories - Dante Sanaei');
                plot([min(xlim()),max(xlim())],[100,100], 'k--')
                
                ylim([0 250])
                coasting_time = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g;
                time_in_space_S1 = [time_in_space_S1 coasting_time]; 
                angles_S1 = [angles_S1 rad2deg(theta)]; 
                if angle == 65 
                    x_optimal_space_S2 = x;
                    y_optimal_space_S2 = y;
                    x1_S2 = x1;
                    y1_S2 = y1;
                    coast_time_OS_S2 = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g;
                    apogee_S2_OS = max(y)
                    x_coast(k)
                end
            end
            if iter == 2
                figure(5)
                plot(x,y); grid on; hold on;
                xlabel('Downrange Position, x (km)');
                ylabel('Altitude, y (km)');
                title('GOLauncher Ballistic Trajectories - Dante Sanaei');
                plot([min(xlim()),max(xlim())],[100,100], 'k--')
                ylim([0 800])
                coasting_time = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g;
                time_in_space_GO = [time_in_space_GO coasting_time];  
                angles_GO = [angles_GO rad2deg(theta)]; 
                if angle == 70 
                    x_optimal_space_GO = x;
                    y_optimal_space_GO = y;
                    x1_GO = x1;
                    y1_GO = y1;
                    coast_time_OS_GO = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g;
                    apogee_GO_OS = max(y)
                    x_coast(k)
                end
            end
        end
    end
    
    figure(6)
    subplot(2,1,1)
    plot(angles_S1, time_in_space_S1); grid on;
    title('SpaceShipTwo Time in Space - Dante Sanaei')
    ylabel('Time (sec)')
    xlabel('Launch Angle (deg)')
    
    subplot(2,1,2)

    plot(angles_GO, time_in_space_GO); grid on    
    title('GOLauncher Time in Space - Dante Sanaei')
    ylabel('Time (sec)')
    xlabel('Launch Angle (deg)')
    
    
    theta = optimal_range_theta;
    Vx0 = V0*cos(theta);
    Vy0 = V0*sin(theta);

    Vx1 = (f*T*cos(theta)) + Vx0;
    Vy1 = (f*sin(theta)-g)*T + Vy0;
    x1 = .5*f*T^2*cos(theta);
    y1 = .5*(f*sin(theta)-g)*T^2 + y_init;
 
    u = sqrt(Vx1^2+Vy1^2);

    x_burn = 0:.1:x1;
    slope = (y1-y_init) / x1;
    y_burn = slope * x_burn + y_init;

    time = 0:.1:final_time;
    x_coast = x1 + Vx1*time ;
    y_coast = y1 + Vy1.*time - .5 * g * time.^2;

    k = find(y_coast >-.01,1, 'last');

    x = [x_burn x_coast(1:k)] / 1000;
    y = [y_burn y_coast(1:k)] / 1000;
    if iter == 1
        figure(7)
        plot(x,y); grid on; hold on
        apogee = max(y)
        xlabel('Downrange Position, x (km)');
        ylabel('Altitude, y (km)');
        title('SpaceShipTwo Optimal Suborbital Trajectory - Dante Sanaei');
        plot(x_optimal_space_S2, y_optimal_space_S2)
        plot([min(xlim()),max(xlim())],[100,100], 'k--'); plot([x1_S2/1000 x1/1000],[y1_S2/1000 y1/1000], 'r*');
        legend('Optimal Range', 'Estimated Optimal Cruise (65 deg)', 'Karman Line', 'Burnout')
       
        coast_time = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g
        max_range_formula = f*T^2 * ( f/g * cot(theta) - .5 * cos(theta))
        max_range_real = x_coast(k)
        coast_time_OS_S2
        rad2deg(optimal_range_theta)
        
        figure(8)
        plot(x,y); grid on; hold on
        xlabel('Downrange Position, x (km)');
        ylabel('Altitude, y (km)');
        S2fg = f/g
        x_iter1 = x1;
        y_iter1 = y1;
        
    end
    if iter == 2
        figure(9)
        plot(x,y); grid on; hold on;
        apogee = max(y)
        xlabel('Downrange Position, x (km)');
        ylabel('Altitude, y (km)');
        title('GOLauncher Optimal Suborbital Trajectory - Dante Sanaei');
        plot(x_optimal_space_GO, y_optimal_space_GO)
        plot([min(xlim()),max(xlim())],[100,100], 'k--'); plot([x1_GO/1000 x1/1000],[y1_GO/1000 y1/1000], 'r*');
        
        coast_time = (sqrt(Vy1^2 + 2 * g * (y1 - 100e3)) + Vy1) / g
        max_range_formula = f*T^2 * ( f/g * cot(theta) - .5 * cos(theta))
        max_range_real = x_coast(k)
        coast_time_OS_GO
        rad2deg(optimal_range_theta)
        
        figure(8)
        plot(x,y); grid on; hold on;
        xlabel('Downrange Position, x (km)');
        ylabel('Altitude, y (km)'); 
        title('Range Optimized Trajectory - Dante Sanaei');
        plot([min(xlim()),max(xlim())],[100,100], 'k--'); plot([x_iter1/1000 x1/1000],[y_iter1/1000 y1/1000], 'r*'); 
        legend('SpaceShipTwo (49.5)', 'GOLauncher (46.87)', 'Karman Line', 'Burnout')
        GOfg = f/g
        
    end
    
    
    
    
    
    
end




%% Functions
function dX_dtau = ascent_odes_tf(tau,X,tf)
    % X(1) = xbar, horizontal component of position
    % X(2) = ybar, vertical component of position
    % X(3) = Vxbar, horizontal component of velocity
    % X(4) = Vybar, vertical component of velocity
    % X(5) = lambda_2_bar, first costate
    % X(6) = lambda_3_bar, second costate
    % X(7) = lambda_4_bar, third costate
    global g Vc h K1 scale F m0 mdot m  
    
    % Adjust mass
    m = m0-abs(mdot)*tau*tf;
    % Simplify long equations
    lam_mag = sqrt(X(6)^2+X(7)^2);
    V_mag = sqrt(X(3)^2+X(4)^2);
    %State and Co-state DE's in terms of d/dt:
    xbardot = X(3)*Vc/h;
    ybardot = X(4)*Vc/h;
    Vxbardot = (F/(m*Vc)) * (-X(6)/lam_mag) - K1*exp(-X(2)*scale) * X(3) * V_mag;
    Vybardot = (F/(m*Vc)) * (-X(7)/lam_mag) - K1*exp(-X(2)*scale) * X(4) * V_mag - (g/Vc);

    if sqrt(X(3)^2+X(4)^2) == 0
        lambda_2_bar = 0;
        lambda_3_bar = 0;
        lambda_4_bar = -X(5)*Vc/h;
    else
        lambda_2_bar = (X(6)*X(3) + X(7)*X(4)) * exp(-X(2)*scale) * (-K1 * scale * V_mag);
        lambda_3_bar = K1 * exp(-X(2) * scale) * (X(6) * (V_mag + (X(3)^2 / V_mag)) + X(7) * (X(3)*X(4)/V_mag));
        lambda_4_bar = (-X(5) * (Vc/h)) + K1 * exp(-X(2) * scale) * (X(7) * (V_mag + (X(4)^2 / V_mag)) + X(6) * (X(3)*X(4)/V_mag));
    end

    dX_dtau = tf*[xbardot; ybardot; Vxbardot; Vybardot;lambda_2_bar; lambda_3_bar; lambda_4_bar];
    return
end

function PSI = ascent_bcs_tf(Y0,Yf,tf)
    global xbar0 ybar0 Vxbar0 Vybar0 ybarf Vxbarf Vybarf Vc F K1 scale g m0 mdot
    Hf = -1;
    mf = m0-abs(mdot)*tf;

    PSI = [Y0(1) - xbar0; % Initial Condition
           Y0(2) - ybar0; % Initial Condition
           Y0(3) - Vxbar0; % Initial Condition
           Y0(4) - Vybar0; % Initial Condition
           Yf(2) - ybarf; % Final Condition
           Yf(3) - Vxbarf; % Final Condition
           Yf(4) - Vybarf; % Final Condition
           (-F/(mf*Vc))*sqrt(Yf(6)^2 + Yf(7)^2) - K1*exp(-scale)*Vc*(Yf(6)*Yf(3))*sqrt(Yf(3)^2) - (Yf(7)*g/Vc) - Hf];
    return
end