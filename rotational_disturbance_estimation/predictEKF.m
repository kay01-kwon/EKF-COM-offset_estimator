function [s_pred, cov_pred] = predictEKF(s_old, cov_old, u, dt, Param)

% Get MOI parameter
J = Param.J;
J_xx = J(1,1);
J_yy = J(2,2);
J_zz = J(3,3);
J_inv = diag([1/J_xx, 1/J_yy, 1/J_zz]);

% Force and moment
f = u(1);
m = [u(2), u(3), u(4)]';

% Prior mean for quaternion
q_prior = [s_old(1);
           s_old(2);
           s_old(3);
           s_old(4)];

qw = q_prior(1);
qx = q_prior(2);
qy = q_prior(3);
qz = q_prior(4);

% Prior mean for angular velocity
w_prior = [s_old(5);
           s_old(6);
           s_old(7)];

wx = w_prior(1);
wy = w_prior(2);
wz = w_prior(3);

% Covert angular velocity to quaternion form
w_quat = [0;
          wx;
          wy;
          wz];

% Derivative of otimes(q,w) w.r.t quaternion
w_prior_skiew = [0,  -wx, -wy, -wz;
                 wx,   0,  wz, -wy;
                 wy, -wz,   0,  wx;
                 wz,  wy, -wx,  0];

% Derivative of otimes(q,w) w.r.t angular velocity
q_prior_skiew = [-qx, -qy, -qz;
                  qw, -qz,  qy;
                  qz,  qw, -qx;
                 -qy,  qx,  qw];

% q_{k+1} = q_{k} + 0.5*otimes(q_{k},w_{k})
% dq_{k+1}/dq_{k} = I + 0.5*d otimes(q_{k},w_{k})/dq_{k}
% dq_{k+1}/dw_{k} = 0.5*d otimes(q_{k},w_{k})/dw_{k}
% dq_{k+1}/dtau_{k} = 0

dq_dq = eye(4) + 0.5*w_prior_skiew*dt;
dq_dw = 0.5*q_prior_skiew*dt;
dq_dtau = zeros(4,3);

dw_dq = zeros(3,4);
dw_dw = dt*J_inv*[0,              (J_zz-J_yy)*wz, (J_zz-J_yy)*wy;
                (J_xx-J_zz)*wz,  0,              (J_xx-J_zz)*wx;
                (J_yy-J_xx)*wy,  (J_yy-J_xx)*wx, 0];
dw_dw = eye(3) + dw_dw;

dw_dtau = dt*J_inv;


dtau_dq = zeros(3,4);
dtau_dw = zeros(3,3);
dtau_dtau = eye(3);


A = [dq_dq, dq_dw, dq_dtau;
     dw_dq, dw_dw, dw_dtau;
     dtau_dq, dtau_dw, dtau_dtau];

tau = [s_old(8), s_old(9), s_old(10)]';


q_pred = q_prior + 0.5*otimes(q_prior, w_quat)*dt;
q_pred = quat_normalize(q_pred);
w_pred = w_prior + J_inv*(m - cross(w_prior,J*w_prior) + tau)*dt;

s_pred = [q_pred;w_pred;tau];

Rt = 0.01^2*eye(10);

cov_pred = A*cov_old*A' + Rt;

end