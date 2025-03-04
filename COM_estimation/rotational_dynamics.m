function sdot = rotational_dynamics(s, u, d, Param)

% s(1) ~ s(4) : qw, qx, qy, qz
% s(5) ~ s(7) : wx, wy, wz

% State variables
% Quaternion and angular velocity
q = [s(1) s(2) s(3) s(4)]';
w = [s(5) s(6) s(7)]';

% Angular velocity in the quaternion form
w_quat = [0;w];

% Inertial Parameter
% Moment of inertia
J = Param.J;
J_xx = J(1,1);
J_yy = J(2,2);
J_zz = J(3,3);
J_inv = diag([1/J_xx, 1/J_yy, 1/J_zz]);

% COM offset
rx = Param.rx;
ry = Param.ry;

% Force and moment
f = u(1);
m = [u(2);u(3);u(4)];

% The effect of COM misalignment
r_cross_f = [ry*f, -rx*f, 0]';

dqdt = 0.5*otimes(q,w_quat);
% Rotational dynamics
dwdt = J_inv*(m - cross(w,J*w) - r_cross_f + d);


sdot = [dqdt;dwdt];

end