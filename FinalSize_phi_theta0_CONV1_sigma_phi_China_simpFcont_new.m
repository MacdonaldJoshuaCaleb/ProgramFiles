load('paramfitDeaths.mat','paramfitDeaths')

ChinaData = readtable('ChinaDataReshaped2.csv');
ChinaData = table2array(ChinaData(:,2:end));
Deaths = readtable('ChinaDeathsReshaped.csv');
Deaths = table2array(Deaths(:,2:end));
MaxDeaths = readtable('MaxDeaths.csv');
MaxDeaths = table2array(MaxDeaths);
MD = sum(MaxDeaths);

for j = 1:90
    CDt(j) = sum(Deaths(j,:));
end

CDt = (CDt./max(CDt))*MD;
CDt = round(CDt);
 CCt = [0;ChinaData(:,end-2)];
 CDt = CDt(1:length(CCt)-1+15);
D = diff(CDt);

C = diff(CCt);

% figure 
% plot(1:70,C)



temp1 = sum(C(3:11));

Slope1 = (C(12)-C(2))./(12-2);
for j = 3:11
    C(j) = C(j-1)+Slope1;
end

Missing1 = temp1 - sum(C(3:11));

temp2 = sum(C(34:36));

Slope2 = (C(33)-C(37))./(33-37);

for j = 34:36
    C(j) = C(j-1)+Slope2;
end

Missing2 = temp2 - sum(C(34:36));

Missing = Missing1 + Missing2;

Cnew = C + Missing*C./CCt(end);
Cnew = round(Cnew);

CCtnew(1) = CCt(1);
for j = 2:length(CCt)
    CCtnew(j) = CCtnew(j-1) + Cnew(j-1);
end
CCtnew = (CCtnew./CCtnew(end)).*(CCt(end));
CCtnew = round(CCtnew);
Cnew = diff(CCtnew);
% figure
% plot(0:70,CCt,'b',0:70,CCtnew,'r')
% 
% figure 
% plot(1:70,C,'b*',1:70,Cnew,'r*')

C = Cnew;
tc = 1:length(C);

td = 0:length(CDt)-1;
tdd = 1:length(D);
[mx,im]=max(C);
ids = 1:length(C);

CC = CCt;
CC = CC(CC > 0);
for i=1:length(C)
    CC(i+1)=CC(i)+C(i);
end



T=4.64;
T_e=2.71;
tau=3;

alphac=1/14;
p = .06;
Mc0 = 0.00172837295386061;
ProvincePops = readtable('ChinaProvincePops.csv');
ProvincePops = table2array(ProvincePops);
N = sum(ProvincePops); % China pop

QCt = ChinaData(:,end);
M=ChinaData(:,end-1);
M = M';
Slope3 = (M(16)-M(3))./(16-3);
for k = 4:15
    M(k) = M(k-1)+Slope3;
end


tq = find(M > 0);

deathdist = @(tv,ag,bg) gamcdf(tv,ag,bg);

param=paramfitDeaths;
beta = param(2)/T;
    sigma = param(3);
    phi = param(4);
    rho = param(5);
    xi = param(6);
    ag = param(7);
    bg = param(8);
    X0 = IC(param,QCt,C,C,M,tq,N,1);

    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*X(1)*((beta/N)*X(3)),(beta/N)*X(1)*X(3)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td,X0);
   %Cdd= zeros(length(tc),1);
   %Cdd(1)=C(1);
   %Cdd(2:end)=Sol1(tc(2:end),7)+Sol1(tc(2:end),8)-(Sol1(tc(1:end-1),7)+Sol1(tc(1:end-1),8));
   tmp = diff(Sol1(:,7))+diff(Sol1(:,8));
   Cdd = tmp(ids);
   Cdd(1) = C(1);
   Ddd = zeros(length(td),1);
   Ddd(1) = CDt(1);
   infs = zeros(length(td),1);
   infs(1) = Sol1(1,end);
   infs(2:end) = Sol1(2:end,end)-Sol1(1:end-1,end);
   for s = 2:length(td)
   mdd(s-1)=(deathdist((2:s),ag,bg)-deathdist((1:s-1),ag,bg))*(infs(s-1:-1:1));
            Ddd(s) = Ddd(s-1)+mdd(s-1);
   end
   Ddd = xi*Ddd;
   Dds = diff(Ddd);
   Dds(1) = D(1);
   %Ddd = Ddd(tdd);
   CumSq=(Sol1(tq,9));
   CumSq(1) = M(tq(1));
   Cds = Cdd(ids);
%    Cumcc=Sol1(tds,7)+Sol1(tds,8);
%    Cumcc(im-1:end)=Cumcc(im-1:end)+C2ex;paramfit(3) + ((1-p)/p)*paramfit(4)*paramfit(5)
  
   z=[Cds',CumSq',Dds'];

  
 Squ=Sol1(:,10);

figure(1)
plot(td,Squ/N*100,'b-')
hold on

figure(2)
plot(ids,Cds,'r-')
hold on
param(4)=0;
 phi = param(4);
    rho = param(5);
    xi = param(6);
    ag = param(7);
    bg = param(8);
    X0 = IC(param,QCt,C,C,M,tq,N,1);

    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*X(1)*((beta/N)*X(3)),(beta/N)*X(1)*X(3)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td,X0);
   %Cdd= zeros(length(tc),1);
   %Cdd(1)=C(1);
   %Cdd(2:end)=Sol1(tc(2:end),7)+Sol1(tc(2:end),8)-(Sol1(tc(1:end-1),7)+Sol1(tc(1:end-1),8));
   tmp = diff(Sol1(:,7))+diff(Sol1(:,8));
   Cdd = tmp(ids);
   Cdd(1) = C(1);
   Ddd = zeros(length(td),1);
   Ddd(1) = CDt(1);
   infs = zeros(length(td),1);
   infs(1) = Sol1(1,end);
   infs(2:end) = Sol1(2:end,end)-Sol1(1:end-1,end);
   for s = 2:length(td)
   mdd(s-1)=(deathdist((2:s),ag,bg)-deathdist((1:s-1),ag,bg))*(infs(s-1:-1:1));
            Ddd(s) = Ddd(s-1)+mdd(s-1);
   end
   Ddd = xi*Ddd;
   Dds = diff(Ddd);
   Dds(1) = D(1);
   %Ddd = Ddd(tdd);
   CumSq=(Sol1(tq,9));
   CumSq(1) = M(tq(1));
   Cds = Cdd(ids);
%    Cumcc=Sol1(tds,7)+Sol1(tds,8);
%    Cumcc(im-1:end)=Cumcc(im-1:end)+C2ex;paramfit(3) + ((1-p)/p)*paramfit(4)*paramfit(5)
  
   z=[Cds',CumSq',Dds'];

  
 Squ=Sol1(:,10);

figure(1)
plot(td,Squ/N*100,'b--')
hold on

figure(2)
plot(tdd,Dds,'r--')
hold on

param(4)=paramfitDeaths(4);
sigma=param(3)*.1;
param(3)=sigma;
 phi = param(4);
    rho = param(5);
    xi = param(6);
    ag = param(7);
    bg = param(8);
    X0 = IC(param,QCt,C,C,M,tq,N,1);

    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*X(1)*((beta/N)*X(3)),(beta/N)*X(1)*X(3)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td,X0);
   %Cdd= zeros(length(tc),1);
   %Cdd(1)=C(1);
   %Cdd(2:end)=Sol1(tc(2:end),7)+Sol1(tc(2:end),8)-(Sol1(tc(1:end-1),7)+Sol1(tc(1:end-1),8));
   tmp = diff(Sol1(:,7))+diff(Sol1(:,8));
   Cdd = tmp(ids);
   Cdd(1) = C(1);
   Ddd = zeros(length(td),1);
   Ddd(1) = CDt(1);
   infs = zeros(length(td),1);
   infs(1) = Sol1(1,end);
   infs(2:end) = Sol1(2:end,end)-Sol1(1:end-1,end);
   for s = 2:length(td)
   mdd(s-1)=(deathdist((2:s),ag,bg)-deathdist((1:s-1),ag,bg))*(infs(s-1:-1:1));
            Ddd(s) = Ddd(s-1)+mdd(s-1);
   end
   Ddd = xi*Ddd;
   Dds = diff(Ddd);
   Dds(1) = D(1);
   %Ddd = Ddd(tdd);
   CumSq=(Sol1(tq,9));
   CumSq(1) = M(tq(1));
   Cds = Cdd(ids);
%    Cumcc=Sol1(tds,7)+Sol1(tds,8);
%    Cumcc(im-1:end)=Cumcc(im-1:end)+C2ex;paramfit(3) + ((1-p)/p)*paramfit(4)*paramfit(5)
  
   z=[Cds',CumSq',Dds'];

 
 Squ=Sol1(:,10);

figure(1)
plot(td,Squ/N*100,'b*')
hold off

figure(2)
plot(tdd,Dds,'r*')
hold off
%  Rfit=beta*T
% %T=5
% Refit=(1-phi*rho)*Rfit
% 
%        % T_e=q(3);
%       
% 
% % E0=beta*(1-phi)*T*I0/(1/tau*T+1)+betaM*(1-phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1));
% 
%        % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) phi*C(1) Mc0];
% %  A=[-1/tau,beta*S0*(1-phi);1/tau,-1/T]
% %  FS=expm(A)
% %  FSI=inv(FS)
% %  x00=(FSI^40)*[E0;I0]
% %  
% %  E00=x00(1)
% %  I00=x00(2)
%  %Rmfit=betaf*T_e*thetafit*S0
%  Refit=(1-phi*rho)*Rfit
%        % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
% %  E0=beta*(1-phi)*T*I0/(1/tau*T+1)+betaM*(1-phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1));
% %         Ec0=beta*(phi)*T*I0/(1/tau*T+1)+betaM*(phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1));
% %         
% %     X0=[S0 Sc0 E0 Ec0 I0 Ic0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0]
% %      X0=[S0 Sc0 E0 Ec0 I0 0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0 0];
% %         
%  Rmfit=0;
% %  Refit=(1-phif)*Rfit+phif*Rmfit
%  R=Rfit
% % Rq=betaM*T_e*N
% theta=Rmfit/Rfit
%  S0=N;
% % I0=I0/N;
% 
% %     function [dxdt]=SEIR(t,X)
% %      dxdt(1)=-beta*X(1)*(X(3));
% %      dxdt(2)= beta*X(1)*(X(3))- sigma*X(2);
% %      dxdt(3)= sigma*X(2)-1/T*X(3);
% %      dxdt=dxdt';
% %     end
% 
% Re=Refit;
% I0=X0(3);
% E0=X0(2);
% %     fun=@(x)  log(x)-(Re*(x-1)-(1+psi)*beta*T*(I0+E0));
% %      Xs=bisection(fun,0,1);
% %      
% %      FS1=CC(1)+1/(1+psi)*N*(1-Xs)
% %      
% %      fun=@(x)  log(x)-(R*(x-1)-(1+sigma)*beta*T*(I0+E0));
% %      Xs=bisection(fun,0,1);
% %      
% %      FS1p0=CC(1)+1/(1+sigma)*N*(1-Xs)
% %      
% %      fsred1=1-FS1/FS1p0
% %      
% %      PS1=N/((1+psi)*R)*(log(1/Re)+Re-1)+(I0+E0)
% %      
% %      PS1p0=N/((1+sigma)*R)*(log(1/(R))+(R)-1)+(I0+E0)
% %      
% %      psred1=1-PS1/PS1p0
% 
% 
% load('parameters_simple_lesst.mat');
% filename = 'data-preview_cn_apr.xlsx';
% xlRange = 'E34:BJ37';
% dataArray = xlsread(filename,xlRange);
% CC=dataArray(1,:);
% QC=dataArray(2,:);
% C=diff(CC);
% Cd=diff(CC);
% Nd=length(Cd);
% td=1:length(CC)
%  N=1e+08;
%  tau=3;
%  T=4.64;
% %T=5
% T_e=2.71;
% T_c=2.71;
%  psi=paramfit(3)
%       % beta=paramfit(4)
%       %T=paramfit(8)
%         p=paramfit(5)
%         phi=paramfit(4)
%        
%       rho=1;
%        % T_e=q(3);
%          sigmaf=psi-(1-p)/p*phi*rho
%         I0=paramfit(1)
%         phif=phi;
%         %Sc0=param(2);
%         S0=N
% Mc0=paramfit(2)
% % beta=betaf;
%  beta=6/N/T;
%  Td=.1609*C;
%  %T_e=Tefit;
%  %T=Tfit;
%   Rfit=beta*T*S0
% 
% % E0=beta*(1-phi)*T*I0/(1/tau*T+1)+betaM*(1-phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1));
%        % Ec0=beta*(phi)*T*I0/(1/tau*T+1)+betaM*(phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1));
%  E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N
%         Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
%         Ic0=phi*rho*I0;
%         X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0 0];
%        % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) phi*C(1) Mc0];
% %  A=[-1/tau,beta*S0*(1-phi);1/tau,-1/T]
% %  FS=expm(A)
% %  FSI=inv(FS)
% %  x00=(FSI^40)*[E0;I0]
% %  
% %  E00=x00(1)
% %  I00=x00(2)
%  %Rmfit=betaf*T_e*thetafit*S0
%  Refit=(1-phi*rho)*Rfit
%        % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
% %  E0=beta*(1-phi)*T*I0/(1/tau*T+1)+betaM*(1-phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1));
% %         Ec0=beta*(phi)*T*I0/(1/tau*T+1)+betaM*(phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1));
% %         
% %     X0=[S0 Sc0 E0 Ec0 I0 Ic0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0]
% %      X0=[S0 Sc0 E0 Ec0 I0 0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0 0];
% %         
%  Rmfit=0;
% %  Refit=(1-phif)*Rfit+phif*Rmfit
%  R=Rfit
% % Rq=betaM*T_e*N
% theta=Rmfit/Rfit
%  S0=N;
% % I0=I0/N;
% 
% %     function [dxdt]=SEIR(t,X)
% %      dxdt(1)=-beta*X(1)*(X(3));
% %      dxdt(2)= beta*X(1)*(X(3))- sigma*X(2);
% %      dxdt(3)= sigma*X(2)-1/T*X(3);
% %      dxdt=dxdt';
% %     end
% 
% % tn=0:100;
% % tnt=1:200;
% % %     function [dxdt]=SEIR(t,X)
% % %      dxdt(1)=-beta*X(1)*(X(3));
% % %      dxdt(2)= beta*X(1)*(X(3))- sigma*X(2);
% % %      dxdt(3)= sigma*X(2)-1/T*X(3);
% % %      dxdt=dxdt';
% % %     end
% % phi1=phi;
% % %phi1=phi;
% %  [tn1,Soln] = ode45(@(t,X)EodeSimpUN(t,X,beta,psi,phi1,p,rho,T),td,X0);
% %  phifit=phi1;
% %  
% %  Csfit=Soln(2:end,7)/rho+Soln(2:end,8)-(Soln(1:end-1,7)/rho+Soln(1:end-1,8));
% % 
% % figure(1)
% % plot(td(2:end)+10,Csfit/N*100,td+10,Soln(:,10)/N*100)
% %  
% % %  peak=max(Csfit)
% % %  
% % % phi_c=min(1,1/(1-theta)*(1-1/R))
% % % %phivec=linspace(0,phi_c,1000);
% % % phi2=0;
% % % sigma2=sigmaf;
% % % phi=phi2;
% % %  E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N
% % %         Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
% % %         Ic0=phi*rho*I0;
% % %         X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
% % % psi2=(1-p)/p*phi2*rho+sigma2
% % %  [tn2,Soln] = ode45(@(t,X)EodeSimpUN(t,X,beta,psi2,phi2,p,rho,T),tnt,X0);
% % %  Csfit=Soln(2:end,7)/rho+Soln(2:end,8)-(Soln(1:end-1,7)/rho+Soln(1:end-1,8));
% % %  
% % %  peak0=max(Csfit)
% % %  
% % %  pred=1-peak/peak0
% % 
% 
% 
% 
