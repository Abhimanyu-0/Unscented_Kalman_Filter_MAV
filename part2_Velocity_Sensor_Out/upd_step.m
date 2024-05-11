function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
 
    % Defined transformations so can be used later
    T_camera_to_robot = [0.707, -0.707, 0, -0.04;
                        -0.707, -0.707, 0, 0;
                         0 , 0,  -1, -0.03;
                          0,0,0,1];

    R_camera_to_robot = [0.707 -0.707 0;
                        -0.707 -0.707 0; 
                          0 0 -1]; 
    R_robot_to_camera = R_camera_to_robot';

    % Defining the parameters for Uscented Transform
    alpha = 0.001;
    k = 1;
    beta = 2;
    n = 15;
    lambda_dash = alpha^2*(n+k)-n;
    covarEstRoot_form = chol(covarEst,"lower");
   
    R = 0.0001*eye(3);

    sigma = [];

    % Calculating sigma points
    for num = 0: 2*n
        
        if (num == 0)

            x = uEst;

        elseif(num > 0 && num <= n)

            x = uEst + (sqrt(n + lambda_dash)*(covarEstRoot_form(:,num)));

        elseif(num > n)
            x = uEst - (sqrt(n + lambda_dash)*(covarEstRoot_form(:,num-n)));

        end
        sigma = [sigma x];
        
    end

 
    ZT = [];
  

    for len = 1: length(sigma)

       angles = [sigma(6, len) sigma(5,len) sigma(4, len)];
       R_body_to_world=(rotz(angles(1))*roty(angles(2))*rotx(angles(3)))';
       skew_matrix = [0 , -T_camera_to_robot(3,4), T_camera_to_robot(2,4);
                   T_camera_to_robot(3,4), 0 , -T_camera_to_robot(1,4);
                   -T_camera_to_robot(2,4),0 , T_camera_to_robot(2,4) ];

       % converted the data to camera frame
       z = R_camera_to_robot*R_body_to_world*[sigma(7:9,len)]-R_camera_to_robot*skew_matrix*R_camera_to_robot*z_t(4:6,1);

           ZT = [ZT z];
    
        end
     p = [];
     for c = 1: length(ZT)
    
        if(c == 1)
            Wm = lambda_dash/(n + lambda_dash);
            p = Wm*ZT(:,c);
           
        else
            Wm = 1/(2*(n + lambda_dash));
            p = p + (Wm*ZT(:,c));
       
        end
    end
    zut = p;

    for c = 1: length(ZT)
    
        if(c == 1)
            Wi = lambda_dash/(n + lambda_dash)+ (1 - alpha^2 + beta);
            C = Wi*(sigma(:,c)-uEst)*transpose(ZT(:,c)-zut);
            S = Wi*(ZT(:,c)-zut)*transpose(ZT(:,c)-zut) + R;
        else
            Wi = 1/(2*(n + lambda_dash));
            C = C + Wi*(sigma(:,c)-uEst)*transpose(ZT(:,c)-zut);
            S = S + Wi*(ZT(:,c)-zut)*transpose(ZT(:,c)-zut) + R;
        end
    
    end
    Ct = C;
    St = S ;
    Kt = Ct*inv(St);
    
    uCurr = uEst + Kt*(z_t(1:3)- zut);
    covar_curr = covarEst - Kt*St*transpose(Kt);
    



end

