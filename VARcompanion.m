function Xi = VARcompanion(A,p,const)

n=size(A,1); %number of observations
%construct companion matrix
Xi=[A(:,const+1:end);eye(n*(p-1)) zeros(n*(p-1),n)];

end

