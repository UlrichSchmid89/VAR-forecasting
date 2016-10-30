clear 
close all
%% Example used in class FS16
%load & fit data 

%load bullionist_controversy



%temp=Year<1824; %temp=use only once, 1 if <1824, everything else 0
%CPI=CPI(temp,:); %CPI=oldCPI for temp=1
%Industrial_Production=Industrial_Production(temp,:);
%Year=Year(temp,:);

clear temp

%% Data used for assignment FS16
% For US data, implement structural breaks: 1-120, 121-173, 174-252

load('USA2.mat');
Industrial_Production = USA2(:,1);
CPI = USA2(:,2);

%%
p=3; %amount of lags, use AIC/BIC to estimate optimal laglength, default was 6. changed to 3

%take logs
lCPI=log(CPI);
lIndustrial_Production=log(Industrial_Production);

dCPI=lCPI(2:end)-lCPI(1:end-1); %create inflation rate, care: we loose one observation with this, end=auotomatically operate from 1 until last
dIP=lIndustrial_Production(2:end)-... %... means, command will continue on next line
    lIndustrial_Production(1:end-1);

%figure(1) 
%plot(Year(2:end),[dCPI dIP]); %2:end because we lost one observation

data1=[dIP';dCPI']; %growth rate of .. production, change such that variable in rows and observations columns + order(p lags means loose p observations)
Y=data1(:,p+1:end);
const=1;
n=size(data1,1);

Z=zmat(data1,p,const);

[A,Sigma]=VAROLS(data1,p,const); 

Xi=VARcompanion(A,p,const);

J=jmat(n,p);

B1=J*inv(eye(n*p)-Xi)*J';   %reduced form long-run multiplier
temp=B1*Sigma*B1;
C1=chol(temp,'lower');   % structural long-run multiplier

S=blanchard_quah(Xi,J,n,p,Sigma);    %impact matrix
epsilon=eye(n);
h=30;
C=zeros(n,n,h+1);

C_diff=impulse_response(Xi,J,S,epsilon,n,h);

C_level=cumsum(C_diff,3);


variables={'y';'p'};
epsilon={'\epsilon_S';'\epsilon_D'};
screen_size=get(0,'Screensize');
f1=figure(1);
set(f1,'Position',[0 0 screen_size(3) screen_size(4)]);
hh=1;
for jj=1:n
    for kk=1:n
        subplot(n,n,hh)
        plot((0:h)',reshape(C_level(jj,kk,:),h+1,1))
        box off
        set(gca,'TickDir','out')
        xlabel('Horizon')
        ylabel('Response')
        title([epsilon{kk} '\rightarrow' variables{jj}])
        line([0 h],[C1(jj,kk) C1(jj,kk)],'LineStyle','--',...
            'Color','k')        
        hh=hh+1;
    end
end


theta=zeros(n,n,h);
theta(:,:,1)=(C_level(:,:,1).^2)./...
    repmat(sum(C_level(:,:,1).^2,2),1,n);
   
for jj=1:h
    theta(:,:,jj+1)=sum(C_level(:,:,1:jj+1).^2,3)./... 
    repmat(sum(sum(C_level(:,:,1:jj+1).^2,3),2),1,n); % 2x2 matrix in our case
end

screen_size=get(0,'Screensize');

f2=figure(2);

set(f2,'Position',[0 0 screen_size(3) screen_size(4)]);

theta1=C1.^2./repmat(sum(C1.^2,2),1,2);

hh=1;
for jj=1:n
    for kk=1:n
        subplot(n,n,hh)
        plot((0:h)',reshape(theta(jj,kk,:),h+1,1))
        set(gca,'TickDir','out')
        box off
        xlabel('Horizon')
        ylabel('Contribution')
        title([epsilon{kk} '\rightarrow' variables{jj}])
        line([0 h], [theta1(jj,kk) theta1(jj,kk)], 'LineStyle','--','Color','k')
        hh=hh+1;
    end
end


filename1='impulse_response_GD.pdf';
%filename1=project_paths('OUT_FIGURES','impulse_response.pdf');
set(f1,'PaperOrientation','landscape');
set(f1,'PaperUnits','normalized','PaperPosition',[0 0 1 1])
print(f1,'-dpdf', filename1)

filename2='fevd_GD.pdf';
%filename2=project_paths('OUT_FIGURES','fevd.pdf');
set(f2,'PaperOrientation','landscape');
set(f2,'PaperUnits','normalized','PaperPosition',[0 0 1 1])
print(f2,'-dpdf', filename2)