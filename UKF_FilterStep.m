function [x_est, C_est,l] = UKF_FilterStep(x, C, measurement, measurementNoiseMean, ...
    measurementNoiseCovariance,l)

measurement=[measurement;1];

x_ukf = [x; measurementNoiseMean];
C_ukf = blkdiag(C, measurementNoiseCovariance);
n = size(x_ukf, 1);
n_state = size(x, 1);
zSigmaPredict=zeros(3,19);

alpha = 1;
beta = 2;
kappa= 0; 
lamda = alpha^2 * (n + kappa) - n;

WM(1) = lamda / (n + lamda);
WM(2 : 2 * n + 1) = 1 / (2 * (n + lamda));

WC(1) = (lamda / (n + lamda)) + (1 - alpha^2 + beta);
WC(2 : 2 * n + 1) = 1 / (2 * (n + lamda));

A = sqrt(n + lamda) * chol(C_ukf)';

xSigma = [zeros(size(x_ukf)) -A A];
xSigma = xSigma + repmat(x_ukf, 1, size(xSigma, 2));
xx = xSigma(1:n_state,:);
noise = xSigma(n_state + 1:n, :);

z = 0;
Sk = 0;
Ck = 0;
B = [l(1)  0  0
     0  l(2)  0
     0   0    0];
  
for i = 1:2 * n+1
    s = noise(1, i);
    v = [noise(2:3, i);0];
    cx = xx(:, i);
    m = [cx(3) cx(5) 0]';
    theta = atan2(measurement(2) - m(2), measurement(1) - m(1));
    the = atan2(cx(2),cx(1));
    xxx = sqrt(l(1)^2 * sin(theta-the) * sin(theta-the) + l(2)^2 * cos(theta-the) * cos(theta-the));
    H = [l(2) * cos(theta-the)/xxx;l(1) * sin(theta-the)/xxx;0];
    F = [cx(1) -cx(2) 0
         cx(2)  cx(1) 0
          0      0    0];    

    zSigmaPredict(:,i) = m + s.*F*B*H + v + [0 0 cx(1)*cx(1)+cx(2)*cx(2)]';% can also be rewritten as a pseudo-measurement
end

for i = 1 : size(zSigmaPredict, 2)
    z = z + ( zSigmaPredict(:,i) * WM(i) );
end

for i = 1 : size(zSigmaPredict, 2)
    Sk = Sk + WC(i) * ( (zSigmaPredict(:,i) - z ) * ( zSigmaPredict(:,i) - z )');
    Ck = Ck + WC(i) * ( (xSigma(1:size(x, 1),i) - x ) * ( zSigmaPredict(:,i) - z )');
end

K = Ck / Sk;
x_est = x + K * (measurement - z);
C_est = C - K * (Sk) * K';
end