function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
    
    %Defining parameters for computing sigma points
 
    alpha = 0.001;
    k = 1;
    beta = 2;
    n = 15;
    n_noise = 12;
    n_prime = n + n_noise;
    lambda_dash = (alpha^2)*(n_prime+k) -n_prime;

    u_augmented = [uPrev; zeros(12,1)];   % 27 X 1
    pAug = [covarPrev zeros(15,12); zeros(12,15) eye(12)*0.01]; 
    pAug_root = chol(pAug,"lower");

    state_aug = [];
    process_model = [];
    xt = [];

    for a = 0: 2*n_prime
        
        if (a == 0)

            x = u_augmented;

        elseif(a > 0 && a <= 27)

            x = u_augmented + (sqrt(n_prime + lambda_dash)*(pAug_root(:,a)));

        elseif(a > 27)
            x = u_augmented - (sqrt(n_prime + lambda_dash)*(pAug_root(:,a-n_prime)));

        end
        state_aug = [state_aug x];
        
    end
    for i = 1:55
    
    
        % Defined for ZYX rotation
    G = [       -sin(state_aug(5,i))                    0        1;
          cos(state_aug(5,i))*sin(state_aug(4,i))    cos(state_aug(4,i))   0;
          cos(state_aug(5,i))*cos(state_aug(4,i))   -sin(state_aug(4,i))   0];

    % Rotation Matrix
    
    R = [cos(state_aug(5,i))*cos(state_aug(6,i)),      cos(state_aug(6,i))*sin(state_aug(4,i))*sin(state_aug(5,i)) - cos(state_aug(4,i))*sin(state_aug(6,i)),      sin(state_aug(4,i))*sin(state_aug(6,i)) + cos(state_aug(4,i))*cos(state_aug(6,i))*sin(state_aug(5,i));
         cos(state_aug(5,i))*sin(state_aug(6,i)),      cos(state_aug(4,i))*cos(state_aug(6,i)) + sin(state_aug(4,i))*sin(state_aug(5,i))*sin(state_aug(6,i)),      cos(state_aug(4,i))*sin(state_aug(5,i))*sin(state_aug(6,i)) - cos(state_aug(6,i))*sin(state_aug(4,i));
                 -sin(state_aug(5,i)),                                   cos(state_aug(5,i))*sin(state_aug(4,i)),                                                     cos(state_aug(4,i))*cos(state_aug(5,i))                               ];               
     
    % gravity
    g = [0;  0; -9.81]; 
    % angular velocity
    wm = [angVel(1,1); angVel(2,1); angVel(3,1)]; 
    % acceleration
    am = [acc(1,1); acc(2,1); acc(3,1)]; 
    % linear velocity
    x3 = [state_aug(7,i); state_aug(8,i); state_aug(9,i)]; 
    % gyroscope bias
    x4 = [state_aug(10,i); state_aug(11,i); state_aug(12,i)]; 
    % accelerometer bias
    x5 = [state_aug(13,i); state_aug(14,i); state_aug(15,i)]; 
    
    % noise definitions
    ng = [state_aug(16,i); state_aug(17,i); state_aug(18,i)];
    
    na = [state_aug(19,i); state_aug(20,i); state_aug(21,i)];
    
    nbg =[state_aug(22,i); state_aug(23,i); state_aug(24,i)];
    
    nba =[state_aug(25,i); state_aug(26,i); state_aug(27,i)];
    
    % Process Model
    
        Xdot = [x3; 
               inv(G)*R*(wm - x4 -ng);
               g+R*(am-x5-na);
               nbg;
               nba];
    
    process_model = [process_model Xdot];
    
        
        %xdot(:,i)
        Xt = (Xdot*dt) + state_aug(1:15, i);
    
        xt = [xt Xt];
    
    end
    

    
    for j = 1: length(xt)
    
        if(j == 1)
            W_mean = lambda_dash/(n_prime + lambda_dash);
            mean = W_mean*xt(:,j);
        else
            W_mean = 1/(2*(n_prime + lambda_dash));
            mean = mean + (W_mean*xt(:,j));
        end
    end
    uEst = mean;
    
    for j = 1: length(xt)
    
        if(j == 1)
            W_cov = (lambda_dash/(n_prime + lambda_dash)) + (1 - alpha^2 + beta);
            cov = W_cov*(xt(:,j)-uEst)*transpose(xt(:,j)-uEst);
        else
            W_cov =  1/(2*(n_prime + lambda_dash));
            cov = cov + (W_cov*(xt(:,j)-uEst)*transpose(xt(:,j)-uEst));
        end
    
    end
    covarEst = cov;



  
end

