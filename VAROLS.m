function [A,Sigma] = VAROLS(data,p,const)

Y=data(:,p+1:end);
Z=zmat(data,p,const);
A=Y*Z'/(Z*Z'); %OLS formula from Script 2.1
T=size(data,2); %number of observations
n=size(data,1); %number of variables
Sigma=(Y*Y'-(Y*Z'/(Z*Z'))*Z*Y')/(T-n*p-const); %formula 2.5 from Script 
                                               %for Variance Covariance Matrix
end

