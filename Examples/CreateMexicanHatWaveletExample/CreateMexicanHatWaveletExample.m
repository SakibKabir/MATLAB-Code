%% Mexican Hat Wavelet
% Create a Mexican hat wavelet with support on [-5,5]. Use 1,000 sample
% points. Plot the result.

% Copyright 2015 The MathWorks, Inc.

lb = 0;
ub = 512;
N = 2^16;
[psi,xval] = mexihat(lb,ub,N);
figure, plot(xval,psi)
title('Mexican Hat Wavelet');
