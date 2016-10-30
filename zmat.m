function Z = zmat(data,p,const)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n=size(data,1); %number of variables
T=size(data,2); %number of observations/time periods

Z=zeros(n*p,T-p); %predefine Z as a matrix of nec. size of zeros

for jj=1:size(Z,2) %build Z-matrix/fill in values
    temp=data(:,jj:jj+p-1); %take relevant part of data for each part of Z
    temp=fliplr(temp); %reverse column ordering (need vector of that)
    Z(:,jj)=temp(:); %inserts selected part stocked as a column into Z
end
Z=[ones(const,size(Z,2));Z];
end

