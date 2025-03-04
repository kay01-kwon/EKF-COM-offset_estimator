close all
clear all


dt = 0.01;
Tf = 5;
time = 0:dt:Tf;
N = length(time);



s_init = [1, 0, 0, 0, ...
          0, 0, 0]';
s_true_vec = zeros(7,N);

Param.J = 0.02*eye(3);
Param.rx = 0.020;
Param.ry = 0;

u = [1, 0., 0, 0]';
d = [0, 0, 0]';

D = 0.1;
w_freq = 1;


s_prior = [1, 0, 0, 0, ...
           0, 0, 0, ...
           0, 0]';

s_est_vec = zeros(9, N);
s_est_vec(:,1) = s_prior;
Cov_prior = 0.01^2*eye(9);
% Cov_prior(8:9,8:9) = 1*eye(2);
% Cov_prior(11:13, 11:13) = 0.1*eye(3);

K_p = 5*eye(3);
K_d = 1*eye(3);

for i = 1:N-1

    % d = [D*cos(w_freq*time(i))-D*w_freq^2, 0, 0]';
    
    s = rk4_ode(@(t,s) rotational_dynamics(s, u, d, Param), s_init, time(i+1), time(i));
    
    q_nom = quat_normalize(s(1:4));
    s(1:4) = q_nom;
    s_true_vec(:,i+1) = s;
    s_init = s;

    dt = time(i+1) - time(i);

    obs = s;

    [s_pred, Cov_pred] = predictEKF(s_prior, Cov_prior, u, dt, Param);
    [s_posterior, Cov_posterior] = measModelEKF(s_pred, Cov_pred, obs);

    s_est_vec(:,i+1) = s_posterior;

    Cov_prior = Cov_posterior;
    s_prior = s_posterior;


    M = pd_control(s, K_p, K_d);
    u(2:4) = M;


end

% Angular velocity plot
figure(1)
subplot(3,1,1)
plot(time, s_true_vec(5,:))
hold on
plot(time, s_est_vec(5,:),'--', LineWidth=2)

subplot(3,1,2)
plot(time, s_true_vec(6,:))
hold on
plot(time, s_est_vec(6,:),'--', LineWidth=2)

subplot(3,1,3)
plot(time, s_true_vec(7,:))
hold on
plot(time, s_est_vec(7,:),'--', LineWidth=2)


% Quaternion plot
figure(2)
subplot(4,1,1)
plot(time, s_true_vec(1,:))
hold on
plot(time, s_est_vec(1,:),'--', LineWidth=2)

subplot(4,1,2)
plot(time, s_true_vec(2,:))
hold on
plot(time, s_est_vec(2,:),'--', LineWidth=2)

subplot(4,1,3)
plot(time, s_true_vec(3,:))
hold on
plot(time, s_est_vec(3,:),'--', LineWidth=2)

subplot(4,1,4)
plot(time, s_true_vec(4,:))
hold on
plot(time, s_est_vec(4,:),'--', LineWidth=2)

% COM offset estimation
figure(3)
subplot(2,1,1)
plot(time, s_est_vec(8,:))

subplot(2,1,2)
plot(time, s_est_vec(9,:))
