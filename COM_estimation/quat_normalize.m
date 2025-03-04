function q_normalized = quat_normalize(q)
%QUAT_NORMALIZE 이 함수의 요약 설명 위치
%   자세한 설명 위치
qw = q(1);
qx = q(2);
qy = q(3);
qz = q(4);

den = sqrt(qw^2 + qx^2 + qy^2 + qz^2);

q_normalized = 1/den*[qw;qx;qy;qz];

end

