function M = pd_control(state,P_gain, D_gain)
%PD_CONTROL 이 함수의 요약 설명 위치
%   자세한 설명 위치
qw = state(1);
q_vec = state(2:4);

w = state(5:7);

M = - P_gain*sign(qw)*q_vec - D_gain*w;

end

