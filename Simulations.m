% clear
% close all
% clc
% 
% tol = 1e-15;
% 
% A = [0.5,2;0,0.9];
% d = [1;0];
% 
% Lambda = [eye(2);-eye(2)];
% lambda = ones(4,1)*2;
% 
% kMax = 5000;
% 
% Comp = zeros(length(lambda),kMax);
% Comp(:,1) = abs(Lambda*d);
% 
% k = 1;
% Ak = A*d;
% 
% while k<kMax && max(Comp(:,k))>tol
%     Comp(:,k+1) = abs(Lambda*Ak);
%     Ak = A*Ak;
%     k = k+1;
% end
% 
% Cand = max(sum(Comp,2));
% clear Comp
% 
% LambdaC = [Lambda,zeros(size(lambda));
%             zeros(2),[Cand;-1]];
% lambdaC = [lambda;1;0];
% 
% [sLambda,slambda] = MRPIscalarDist(LambdaC,lambdaC,A,d,1000,'Moritz');
% 
% % -------
% % Example 3 in Table 1  Schulze-Darup, Schaich & Cannon:
% 
% A = [.5,0;
%     .5,1];
% B = [1;0];
% K = dlqr(A,B,eye(2),1e-2);
% A = A-B*K;
% D = -eye(2);
% 
% Lambda = [-K;
%           K;
%           eye(2);
%           -eye(2)];
% lambda = [.2;.3;.7;.5;.3;.5];
% 
% G = [speye(2);-speye(2)]/.05;
% 
% P = A'*A;
% V = [1,-1,1,-1;
%     1,1,-1,-1]/20;
% cand = zeros(1,4);
% for i = 1:4
%     cand(i) = V(:,i)'*P*V(:,i);
% end
% R = sqrt(max(cand));
% 
% decayRate = sqrt(1-max(eig(A)));
% lambdaMin = min(svd(A));
% Fact = R/(1-decayRate)*sqrt(1/lambdaMin);
% cand = zeros(1,length(lambda));
% for i = 1:length(lambda);
%     cand(i) = Fact*sqrt(Lambda(i,:)*P*Lambda(i,:)');
% end
% hTilde = max(cand);
% 
% LambdaLifted = [Lambda,zeros(size(Lambda,1),1);
%                 zeros(2),[-1;1/1.5]];
% lambdaLifted = [lambda;0;1];
% 
% [TermA,TermB,~] = MRPIFullDimDist(LambdaLifted,lambdaLifted,G,A,D,200,'MPC book Example');
% 
% Num = 10;
% InvSetA = cell(1,Num);
% InvSetb = cell(1,Num);
% iters = zeros(1,Num);
% 
% alpha = linspace(0,1.5,Num);
% for i = 1:Num
%     [InvSetA{i},InvSetb{i},iters(i)] = MRPIfixedAlpha(Lambda,lambda,A,D,G,200,alpha(i));
% end
% 
% 
% figure(1)
% hold on
% for i = 1:Num
%     plot(Polyhedron(InvSetA{i},InvSetb{i}),'alpha',.2)
% end
% 
% 
% 
% x0 = [.25;-.5];
% x = [x0,zeros(2,50)];
% for i = 1:50
%     x(:,i+1) = A*x(:,i);
% end
% 
% plot(x(1,:),x(2,:),'bx')
% hold off


% -----
% Example Mayne, Seron, Rakovic

A = [1,1;0,1];
B = [.5;1];
D = eye(2);
K = dlqr(A,B,eye(2),1e-2);
A = A-B*K;

Lambda = [K;-K;[0,1]];
lambda = [1;1;2];

G = [speye(2);-speye(2)]*10;

LambdaLifted = [Lambda,zeros(size(Lambda,1),1);
                zeros(2),[-1;approxAlphaStar(Lambda,A,D,G,200,1e-16)]];
lambdaLifted = [lambda;0;1];

% [TermA,TermB,iters] = MRPIFullDimDist(LambdaLifted,lambdaLifted,G,A,D,200,'Mayne Example');

Num = 3;
InvSetA = cell(1,Num);
InvSetb = cell(1,Num);
iters = zeros(1,Num);
V = cell(1,Num);
k = cell(1,Num);
al = [1,1.6849,2.5];
eps = (al(2)-al(3))/(al(1)-al(3));

for i = 1:Num
    [InvSetA{i},InvSetb{i},iters(i)] = MRPIfixedAlpha(Lambda,lambda,A,D,G,200,al(i));
    V{i} = vertexEnumeration(InvSetA{i},InvSetb{i});
    k{i} = convhull(V{i}(:,1),V{i}(:,2));
end

CM = colormap;
figure(1)
fill(V{1}(k{1},1),V{1}(k{1},2),CM(1,:),'FaceAlpha',1)
xlabel('$x_1$')
ylabel('$x_2$')

figure(2)
fill(V{2}(k{2},1),V{2}(k{2},2),CM(1,:),'FaceAlpha',1)
Vinter = eps*repmat(V{1},size(V{3},1),1)+(1-eps)*kron(ones(size(V{1},1),1),V{3});
kinter = convhull(Vinter(:,1),Vinter(:,2));
hold on
plot(Vinter(kinter,1),Vinter(kinter,2),'-.','linewidth',2)
hold off
xlabel('$x_1$')
ylabel('$x_2$')

figure(3)
fill(V{3}(k{3},1),V{3}(k{3},2),CM(1,:),'FaceAlpha',1)
xlabel('$x_1$')
ylabel('$x_2$')