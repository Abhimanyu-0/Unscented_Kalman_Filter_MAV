function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    Ct = [1, 0, 0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
          0,1,0,0,0,0,0,0,0,0, 0,0 ,0 ,0 ,0;
          0,0,1,0,0,0,0,0,0,0, 0,0 ,0 ,0 ,0;
          0,0,0,1,0,0,0,0,0,0, 0,0 ,0 ,0 ,0;
          0,0,0,0,1,0,0,0,0,0, 0,0 ,0 ,0 ,0;
          0,0,0,0,0,1,0,0,0,0, 0,0 ,0 ,0 ,0];

     R = [0.00001, 0, 0, 0, 0, 0;
            0, 0.00001, 0, 0, 0, 0;
            0, 0, 0.00001, 0, 0, 0;
            0, 0, 0, 0.00001, 0, 0;
            0, 0, 0, 0, 0.00001, 0;
            0, 0, 0, 0, 0, 0.00001];

     Kalman = covarEst*(Ct')/((Ct*covarEst*Ct' + R));
     
     uCurr =(uEst + Kalman * (z_t - Ct * uEst));
     covar_curr = covarEst - Kalman * (Ct*covarEst);
    

end

