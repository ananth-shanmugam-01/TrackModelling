function [TrackData] = createTrackModel(fileName,trackName,source)

    TrackData = struct;
    TrackData.Name = trackName;
    TrackData.source = source;

    % Allow to choose calculation type

    % switch type

    % case GPS
    rawData = load(fileName);
    time = rawData.data_lap2.t;
    measVelocity = rawData.data_lap2.Chassis_Speed_mps;
    latData = smoothdata(rawData.data_lap2.GPS_Latitude_deg,"gaussian","SmoothingFactor",0.0025);
    longData = smoothdata(rawData.data_lap2.GPS_Longitude_deg,"gaussian","SmoothingFactor",0.0025);
    altData = rawData.data_lap2.GPS_Altitude_m;
    
    [kt, calcSpeed, dist] = GPScurvature(latData,longData,altData,time);
    
    trackDist = 0:sector_dist:dist(end);
    curvatureSpline = csaps(dist, kt,0.8,trackDist);
        
    figure
    hold on
    plot(dist,kt,'DisplayName','GPS')
    plot(trackDist,curvatureSpline,'DisplayName','Interp')
    hold off

    % case Lateral Acceleration
    
    
    function [kt, calcSpeed, dist] = GPScurvature(latData,longData,altData,time)
    
    origin = [latData(1),longData(1),altData(1)];
    [latCart,longCart,~] = latlon2local(latData,longData,altData,origin); % convert to cartesian coords
    
    xdot = gradient(latCart)./gradient(time);
    ydot = gradient(longCart)./gradient(time);
    xddot = gradient(xdot)./gradient(time);
    yddot = gradient(ydot)./gradient(time);
    kt = (xdot.*yddot - ydot.*xddot)./(( xdot.^2 + ydot.^2).^1.5);
    calcSpeed = sqrt(xdot.^2 + ydot.^2);
    dist = cumtrapz(time,calcSpeed);
    
    end

end
