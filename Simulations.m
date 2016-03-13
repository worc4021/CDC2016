clear
close all
clc

tol = 1e-15;

A = [0.5,2;0,0.9];
d = [1;0];

Lambda = [eye(2);-eye(2)];
lambda = ones(4,1)*2;

kMax = 5000;

Comp = zeros(length(lambda),kMax);
Comp(:,1) = abs(Lambda*d);

k = 1;
Ak = A*d;

while k<kMax && max(Comp(:,k))>tol
    Comp(:,k+1) = abs(Lambda*Ak);
    Ak = A*Ak;
    k = k+1;
end

Cand = max(sum(Comp,2));
clear Comp

LambdaC = [Lambda,zeros(size(lambda));
            zeros(2),[Cand;-1]];
lambdaC = [lambda;1;0];

[sLambda,slambda] = MRPIscalarDist(LambdaC,lambdaC,A,d,1000,'Moritz');


