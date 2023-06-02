function fcn_plot_imu_with_frame(sensor_radius, frame_length, wolrld_T_imu)

    % Viusal geometry of the IMU sensors
    [x,y,z] = sphere;
    x       = x * sensor_radius;
    y       = y * sensor_radius;
    z       = z * sensor_radius;   

    % CS of the IMUs ------------------------------------------------------
    trplot(wolrld_T_imu,'rgb', 'notext','length', frame_length,'thick',2)

    % Sphere representing the IMU -----------------------------------------
    sph = surf(x + wolrld_T_imu(1,4), ...
               y + wolrld_T_imu(2,4), ...
               z + wolrld_T_imu(3,4));      
    set(sph,'FaceColor',[1 0 0], ...
      'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none') 
end