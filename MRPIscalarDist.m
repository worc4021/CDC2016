function [Aout,bout] = MRPIscalarDist(Ain,bin,A,D,iMax,chars)

if nargin<6
    chars = ' ';
end

Lambda = Ain;
lambda = bin;

LambdaNext = [];
lambdaNext = [];

i = 1;
while and(~isContained(Lambda,lambda,LambdaNext,lambdaNext),i<iMax)
    
    if mod(i,10) == 0
        clc
        fprintf([chars,'\n']);
        fprintf('Currently at iteration %d, there are currently %d constraints\n',i,length(lambda))
    end
    
    if i~=1
        Lambda = LambdaNext;
        lambda = lambdaNext;
    end
    
    LambdaNext=[Lambda(:,1:2)*A,abs(Lambda(:,1:2)*D)+Lambda(:,3)];
    lambdaNext=lambda;
    
    [LambdaNext,lambdaNext] = bigReduce([Lambda;LambdaNext],...
        [lambda;lambdaNext]);
    i = i+1;

end

Aout = LambdaNext;
bout = lambdaNext;