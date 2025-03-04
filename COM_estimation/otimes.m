function q_result = otimes(q_left, q_right)
% The multiplication of two quaternions

qw_l = q_left(1);
qx_l = q_left(2);
qy_l = q_left(3);
qz_l = q_left(4);

q_l_matrix = [qw_l, -qx_l, -qy_l, -qz_l;
              qx_l, qw_l, -qz_l, qy_l;
              qy_l, qz_l, qw_l, -qx_l;
              qz_l, -qy_l, qx_l, qw_l];

q_result = q_l_matrix*q_right;
end

