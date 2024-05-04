function trackData = initialiseTrackData(trackCoordinates)

    raw_lat_pos     = smoothdata(trackCoordinates.latitude_deg,"gaussian","SmoothingFactor",0.001);
    raw_long_pos    = smoothdata(trackCoordinates.longitude_deg,"gaussian","SmoothingFactor",0.001);
    raw_alt_pos     = trackCoordinates.altitude_m;   
    
    lat_origin  = raw_lat_pos(1);
    long_origin = raw_long_pos(1);
    alt_origin  = raw_alt_pos(1);
    
    origin = [lat_origin, long_origin, alt_origin];
    
    [x, y, z] = latlon2local(raw_lat_pos, raw_long_pos, raw_alt_pos, origin);
    
    % sLap calculation from absolute distance
    dS                  = sqrt(gradient(x).^2 + gradient(y).^2);
    sLap                = cumtrapz(dS);
    
    % HyperSample x,y, sLap
    sLap_fine           = linspace(0,sLap(end),6000);
    x                   = interp1(sLap,x,sLap_fine,"makima");
    y                   = interp1(sLap,y,sLap_fine,"makima");
    
    % Splines
    x_spline            = spline(sLap_fine, x);
    x_der_spline        = fnder(x_spline,1);
    x_der_der_spline    = fnder(x_spline,2);
    
    y_spline            = spline(sLap_fine, y);
    y_der_spline        = fnder(y_spline,1);
    y_der_der_spline    = fnder(y_spline,2);

    % Spline evaluation appears to be smoother when resample coarser
    sLap_coarse         = linspace(0,sLap_fine(end),1000);
    x_dot               = ppval(x_der_spline,sLap_coarse);
    x_ddot              = ppval(x_der_der_spline, sLap_coarse);
    y_dot               = ppval(y_der_spline,sLap_coarse);
    y_ddot              = ppval(y_der_der_spline, sLap_coarse);
    
    % Curvature calculations
    Curv = (x_dot.* y_ddot - x_ddot.*y_dot)./( x_dot.^2 + y_dot.^2).^1.5; 
   
    % Connect Outputs

    trackData = struct;
    trackData.sLap  = sLap_coarse;
    trackData.x     = interp1(sLap_fine,x,sLap_coarse,"makima");
    trackData.y     = interp1(sLap_fine,y,sLap_coarse,"makima");
    trackData.Curv  = csaps(sLap_coarse,Curv,0.3,sLap_coarse);
    trackData.theta = cumtrapz(trackData.sLap, trackData.Curv);
    
    trackData.heading_origin = atan(y_dot(1)/x_dot(1));
    trackData.theta = trackData.theta + trackData.heading_origin;


end