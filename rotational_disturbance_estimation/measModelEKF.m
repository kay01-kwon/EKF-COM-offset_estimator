function [s_est, cov_est] = measModelEKF(s_pred,cov_pred, obs)


Ht = zeros(7,10);
Ht(1:4,1:4) = eye(4);
Ht(5:7,5:7) = eye(3);

Qt = 0.01^2*eye(7);
sensor_model = s_pred(1:7);
innovation = obs - sensor_model;
Kt = cov_pred*Ht'*inv(Ht*cov_pred*Ht'+Qt);
s_est = s_pred + Kt*innovation;
cov_est = (eye(10) - Kt*Ht)*cov_pred;

end