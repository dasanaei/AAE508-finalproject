clc; clear all; close all;
%% PART THREE:  Range Optimal Missile Launch From Ground and Air

% Drag on the vehicle and varying mass are not calculated. 

%% Range of Ground vs Air Missile Launch 
final_time=1000;
for groundair = 1:2
    g = 9.81;
    mdot = 65;
    F = 250e3;
    Ispm = F/mdot;
    m = 5700/1.35;
    f = (mdot * Ispm)/m;
    T = 60;
    g2f = g/f;
    if groundair == 1 
        y_init = 0;
        V0 = 0;
    end
    if groundair == 2
        y_init = 10e3;
        V0 = 650;
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
    
    figure(1)
    plot(x,y); grid on; hold on;
    xlabel('Downrange Position, x (km)');
    ylabel('Altitude, y (km)');
    title('Missile Range Optimized Trajectory - Dante Sanaei');
    hold on; plot(x1/1000, y1/1000, 'r*')
    legend('Ground Launch', 'Burnout', 'Air Launch')
    max_range_formula = f*T^2 * ( f/g * cot(theta) - .5 * cos(theta))
    max_range_real = x_coast(k)
    max(y)
    
end
    



