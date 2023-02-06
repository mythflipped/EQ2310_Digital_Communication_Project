function b = training_sequence(nr_training_bits)
% b = training_sequence(nr_training_bits)
%
% Generate a training sequence consisting of n bits. Currently, a random
% sequence is used, but a deterministic sequence with better
% autocorrelation properties could be used instead.
%
% Input
%   nr_training_bits = length of training sequence (should be an even number)
%
% Output
%   b = training sequence {0/1}
%
% History:
%   2000-09-28  written /Stefan Parkvall
%   2001-10-21  modified / George Jöngren

b = (rand(1, nr_training_bits) > .5);
