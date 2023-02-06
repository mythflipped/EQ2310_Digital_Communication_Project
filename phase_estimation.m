function phihat = phase_estimation(r, b_train)
% phihat = phase_estimation(r, b_train)
%
% Phase estimator using the training sequence. The phase estimate is
% obtained by minimizing the norm of the difference between the known
% transmitted QPSK-modulated training sequence and the received training
% part. NB! There are other ways of estimating the phase, this is just
% one example.
%
% Input:
%   r       = received baseband signal
%   b_train = the training sequence bits
%
% Output:
%   phihat     = estimated phase
phihat = 0;
qpsk_train=qpsk(b_train);
r_train=r(1:length(qpsk_train));
min=norm(r_train-qpsk_train);
for phi=-pi:0.01:pi
    r_phi=r_train*exp(-1j*phi);
    delta=norm(r_phi-qpsk_train);
    if delta<min
        min=delta;
        phihat=phi;
    end
end
end


