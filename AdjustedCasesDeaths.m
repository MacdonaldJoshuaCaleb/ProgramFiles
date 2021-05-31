function [] = AdjustedCasesDeaths
close all

ChinaData = readtable('ChinaDataReshaped.csv');
ChinaData = table2array(ChinaData(:,2:end));
Deaths = readtable('ChinaDeathsReshaped.csv');
Deaths = table2array(Deaths(:,2:end));
MaxDeaths = readtable('MaxDeaths.csv');
MaxDeaths = table2array(MaxDeaths);
MD = sum(MaxDeaths);

for j = 1:79
    CDt(j) = sum(Deaths(j,:));
end

CDt = (CDt./max(CDt))*MD;
CDt = round(CDt);
CDt = CDt(1:57+11);
D = diff(CDt);
CCt = ChinaData(:,end-2);


M=ChinaData(:,end-1);
M = M(2:end)';

Cd2 = diff(CCt);
if Cd2(1) == 0
   %Cd(1) = 2*Cd(2)-Cd(3);
    Cd2(1) = (Cd2(2))/2;
    Cd2(2)= (Cd2(2))/2;
    if Cd2(2) == 0
            Cd2(1) = Cd2(3)/3;
            Cd2(2) = Cd2(3)/3;
            Cd2(3) = Cd2(3)/3;
    end
end
if Cd2(end) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
    Cd2(end) = (Cd2(end-1))/2;
    Cd2(end-1)= (Cd2(end-1))/2;
    if Cd2(end-1) == 0
            Cd2(end) = Cd2(end-2)/3;
            Cd2(end-1) = Cd2(end-2)/3;
            Cd2(end-2) = Cd2(end-2)/3;
    end
end
if Cd2(1) < 0
    Cd2(1) = -Cd2(1);
end
ind0=find(Cd2==0);
   if isempty(ind0)==0
    Cd2(ind0)=1/3*Cd2(ind0-1)+1/3*Cd2(ind0+1);
    Cd2(ind0-1)=2/3*Cd2(ind0-1);
    Cd2(ind0+1)=2/3*Cd2(ind0+1);
   end
   figure 
plot(1:70,Cnew,'r*')
slope = (Cd2(24)-Cd2(20))./(24-20);
temp = sum(Cd2(21:23));

temp3 = Cd2;
Cd2(21) = Cd2(20)+slope;
Cd2(22) = Cd2(20)+2*slope;
Cd2(23) = Cd2(20)+3*slope;
Missing = temp - sum(Cd2(21:23));
Cd2new = Cd2 + Missing*Cd2./CCt(end);
Cd2new = round(Cd2new);

CCtnew(1) = CCt(1);
for j = 2:length(CCt)
    CCtnew(j) = CCtnew(j-1) + Cd2new(j-1);
end
CCtnew = (CCtnew./CCtnew(end)).*(CCt(end));
CCtnew = round(CCtnew)
Cd2new = diff(CCtnew)
figure
hold on
plot(0:57,CCt,'b')
plot(0:57,CCtnew,'r')
hold off

figure
hold on
plot(1:57,Cd2new,'r*')
plot(1:57,temp3,'b*')
hold off


C = Cd2new;
tc = 1:length(C);
tq = find(M > 0);
td = 1:length(D);

[mx,im]=max(C);
ids = 1:length(C);

CC = CCt;
CC = CC(CC > 0);
for i=1:length(Cd2)
    CC(i+1)=CC(i)+Cd2(i);
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
deathdist = @(tv,ag,bg) gamcdf(tv,ag,bg);

function X0 = IC(param,QC,C1,C2,M,tq,N,its)
    I0 = param(1);
    beta = param(2);
    phi = param(4);
    rho = param(5);
    S0=N;
    
    if its == 1
        C = C1;
        
    end
    if its == 2
        C = C2;
    
    end
    E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1));
    Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*(C(1)/(C1(1)+C2(1)));
    Ic0=phi*rho*I0*(C(1)/(C1(1)+C2(2)));
    X0=[S0 E0 I0 Ec0 Ic0 p*QC(1) C(1) phi*rho*C(1) (1-p)*M(tq(1))*(C(tq(1))/(C1(tq(1))+C2(tq(1))))];
    %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end


function[z] = Solve_sys(param,td)
    beta = param(2);
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
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td,X0);
   Cdd= zeros(length(tc),1);
   Cdd(1)=Cd2(1);
   Cdd(2:end)=Sol1(tc(2:end),7)+Sol1(tc(2:end),8)-(Sol1(tc(1:end-1),7)+Sol1(tc(1:end-1),8));
   Ddd = zeros(length(td),1);
   Ddd(1) = D(1);
   tmp = zeros(length(td),1);
   tmp(1) = 0;
   infs = zeros(length(td),1);
   infs(1) = Cd2(1);
   infs(2:end) = Sol1(td(2:end),7)+Sol1(td(2:end),8)-(Sol1(td(1:end-1),7)+Sol1(td(1:end-1),8));
   mdd(1) = D(1);
   for s = 2:length(td)
       mdd(s-1)=(deathdist((2:s),ag,bg)-deathdist((1:s-1),ag,bg))*(infs(s-1:-1:1)+(1-rho)*infs(s-1:-1:1));
       tmp(s) = tmp(s-1)+mdd(s-1);
   end
   tmp = xi*tmp;
   Ddd(2:end) = tmp(2:end)-tmp(1:end-1);
   CumSq=(Sol1(tq-2,9))*sqwt ;
   CumSq(1) = M(tq(1))*sqwt;
   Cds = Cdd(ids)*sqwt2;
%    Cumcc=Sol1(tds,7)+Sol1(tds,8);
%    Cumcc(im-1:end)=Cumcc(im-1:end)+C2ex;
  
   z=[Cds',CumSq',Ddd'];
      
      
end


I01g = (C(1)/.5)*T;
% fit values from entire country fit with exception of I0
paramguess = [I01g,(1.66014398889413e-08)*(1e8),	1407.12998368299,...
              0.391575667314593,.5,.01,6,3.5];
lb = [I01g*.05, 2/T , 0 , .01 ,.01,.001,0,0];
ub = [I01g*10,6/T, inf , 1  ,1,.05,12,inf,inf];
sqwt2 = .1;
sqwt=max(C(ids))./max(M(tq));
length(Solve_sys(paramguess,td))
length([C(ids),M(tq)*sqwt,D])
[paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C(ids)*sqwt2,M(tq)*sqwt,D],lb,ub);
fits = Solve_sys(paramfit,td);
cfit = fits(1:length(ids))./sqwt2;
qfit = fits(length(ids)+1:length(tq)+length(ids))/sqwt;
dfit =  fits(length(tq)+length(ids)+1:end);

array2table(paramfit)
paramfit(end)*paramfit(end-1)
paramfit(2)*T

figure
hold on
plot(ids,(1-paramfit(5))*cfit+cfit,'b-.','linewidth',2)
plot(ids,cfit,'b',ids,C(ids),'k*','linewidth',2)
title('Daily Cases','fontsize',16,'interpreter','latex')
xlabel('Days Since January 23, 2020','fontsize',16,'interpreter','latex')
set(gca,'FontSize',16)
hold off
legend({'True Daily Case Est.','Reported Daily Case Fit'})
xlim([ids(1),ids(end)])

figure
plot(tq,qfit,'b',tq,M(tq),'k*','linewidth',2)
title('Quarintined','fontsize',16,'interpreter','latex')
xlabel('Days Since January 23, 2020','fontsize',16,'interpreter','latex')
set(gca,'FontSize',16)
xlim([tq(1),tq(end)])

figure
plot(td,dfit,'b',td,D,'k*','linewidth',2)
title('Deaths','fontsize',16,'interpreter','latex')
xlabel('Days Since January 23, 2020','fontsize',16,'interpreter','latex')
set(gca,'FontSize',16)
xlim([td(1),td(end)])
ylim([0 inf])
paramfit(3) + ((1-p)/p)*paramfit(4)*paramfit(5)
end