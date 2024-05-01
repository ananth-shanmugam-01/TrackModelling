clear; clc;

load('C:\Users\admin\Desktop\GitHub\QSS-LTS\dataFiles\230722_Endurance_lap2.mat')

lat_pos     = smoothdata(data_lap2.GPS_Latitude_deg,"gaussian","SmoothingFactor",0.001);
long_pos    = smoothdata(data_lap2.GPS_Longitude_deg,"gaussian","SmoothingFactor",0.001);
alt_pos     = data_lap2.GPS_Altitude_m;   

lat_origin  = lat_pos(1);
long_origin = long_pos(1);
alt_origin  = alt_pos(1);

origin = [lat_origin, long_origin, alt_origin];

[x, y, z] = latlon2local(lat_pos, long_pos, alt_pos, origin);


sLap_vCar        = cumtrapz(gradient(data_lap2.t).*data_lap2.Chassis_Speed_mps);

% sLap calculation from absolute distance
dS = sqrt(gradient(x).^2 + gradient(y).^2);
sLap = cumtrapz(dS);

% Curvature calculation

x_dot   = gradient(x)./gradient(sLap);
x_ddot  = gradient(x_dot);
y_dot   = gradient(y)./gradient(sLap);
y_ddot  = gradient(y_dot);

C = (x_dot.* y_ddot - x_ddot.*y_dot)./( x_dot.^2 + y_dot.^2).^1.5; 

figure(1);
plot(sLap, C)
ylim([-0.3, 0.3])



% Heading angle calculation in [rad]




% figure(1);
% hold on
% plot(sLap_vCar,'DisplayName','vCar-based')
% plot(sLap_GPS,'DisplayName','GPS-based')
% hold off;
% legend
% 
% figure(2);
% nexttile
% scatter(x,y,[],sLap_GPS)
% colorbar