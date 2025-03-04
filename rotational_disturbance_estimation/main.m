close all
clear all


dt = 0.01;
Tf = 10;
time = 0:dt:Tf;
N = length(time);



s_init = [1, 0, 0, 0, ...
          0, 0, 0]';
s_true_vec = zeros(7,N);

Param.J = 0.04*eye(3);
Param.rx = 0.020;
Param.ry = 0;

Param_est.J = 0.02*eye(3);

u = [1, 0., 0, 0]';
d = [0, 0, 0]';
d_true_vec = zeros(3,N);
d_true_vec(:,1) = d;


D = 0.1;
w_freq = 1;


s_prior = [1, 0, 0, 0, ...
           0, 0, 0, ...
           0, 0, 0]';

s_est_vec = zeros(10, N);
s_est_vec(:,1) = s_prior;
Cov_prior = 0.01^2*eye(10);
Cov_prior(8:10,8:10) = 1*eye(3);
var_tau_vec = zeros(3,N);
var_tau_vec(:,1) = [Cov_prior(8,8);
                    Cov_prior(9,9);
                    Cov_prior(10,10)];

std_noise = [0.01^2;
             0.01^2;
             0.01^2;
             0.01^2;
             0.01^2;
             0.01^2;
             0.01^2;];

K_p = 3*eye(3);
K_d = 1*eye(3);

rx = Param.rx;
ry = Param.ry;

for i = 1:N-1

    d = [D*cos(w_freq*time(i)) - D*w_freq^2;
        0;
        0];

    if time(i) >= 3
        d(3) = 1;
    end


    r_cross_f = [ry*u(1);
                -rx*u(1);
                0];

    d_true_vec(:,i+1) = d - r_cross_f;
    
    s = rk4_ode(@(t,s) rotational_dynamics(s, u, d, Param), s_init, time(i+1), time(i));
    
    q_nom = quat_normalize(s(1:4));
    s(1:4) = q_nom;
    s_true_vec(:,i+1) = s;
    s_init = s;

    dt = time(i+1) - time(i);

    obs = s + normrnd(0, std_noise);

    [s_pred, Cov_pred] = predictEKF(s_prior, Cov_prior, u, dt, Param_est);
    [s_posterior, Cov_posterior] = measModelEKF(s_pred, Cov_pred, obs);

    s_est_vec(:,i+1) = s_posterior;

    Cov_prior = Cov_posterior;
    s_prior = s_posterior;

    var_tau = [Cov_posterior(8,8);
               Cov_posterior(9,9);
               Cov_posterior(10,10)];

    var_tau_vec(:,i+1) = var_tau;

    d_est = [s_posterior(8), s_posterior(9), s_posterior(10)]';

    w = s(5:7);


    M = pd_control(s, K_p, K_d) - d_est + cross(w,Param.J*w);
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

% lumped disturbance estimation
figure(3)
subplot(3,1,1)
plot(time, d_true_vec(1,:))
hold on
plot(time, s_est_vec(8,:),'--', LineWidth=2)


subplot(3,1,2)
plot(time, d_true_vec(2,:))
hold on
plot(time, s_est_vec(9,:),'--', LineWidth=2)


subplot(3,1,3)
plot(time, d_true_vec(3,:))
hold on
plot(time, s_est_vec(10,:),'--', LineWidth=2)


% Variance of lumped disturbance
figure(4)
subplot(3,1,1)
plot(time, var_tau_vec(1,:))

subplot(3,1,2)
plot(time, var_tau_vec(2,:))

subplot(3,1,3)
plot(time, var_tau_vec(3,:))
