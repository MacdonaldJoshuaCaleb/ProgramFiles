function []= FitSimpleModelChinaDailyVarT
close all
clear 
warning off;
set(0,'DefaultFigureVisible','on')
% stop/start figures from displaying in matlab, switch 'off' to 'on to display
% import data 
ChinaData = readtable('ChinaDataReshaped.csv');
%Province = readtable('ProvinceNames.csv');
%Province = table2array(Province);
%Province = string(Province);
ProvincePops = readtable('ChinaProvincePops.csv');
ProvincePops = table2array(ProvincePops);
% national total cases
CCt = table2array(ChinaData(:,end-2));
%CData2 = table2array(ChinaData(:,15));
%CData1 = CCt - CData2;
% national total quarintined 
QCt = table2array(ChinaData(:,end));
% national current quarintined
%Mt = table2array(ChinaData(:,end));
%Qtd = Mt(2:end);
% national daily cases
%Ctd = diff(CCt);
%tdt=1:length(CCt);
% storage bin for fit parameteres 
%FitparamsAgg = zeros(2,5);




T=4.64;
T_e=2.71;
tau=3;
rho=1;
sqwt=0.05;
alphac=1/14;
p = .06;
%I0 = 1600;
Mc0 = 0.00172837295386061;

%M0=121*exp(r*(Tc(end)-3))/200
N = sum(ProvincePops); % China pop
% retrieve province name
%Name = Province(its);
% provincial total cases
CC = CCt;
td=2:length(CC);
CC = CC(CC > 0);

QC = QCt;
M=table2array(ChinaData(:,end-1));
Qd = M(2:end);
% provincial daily cases
Cd = diff(CC);
% account for no new reported cases on day after first reported case
if Cd(1) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
    Cd(1) = (Cd(2))/2;
    Cd(2)= (Cd(2))/2;
    if Cd(2) == 0
            Cd(1) = Cd(3)/3;
            Cd(2) = Cd(3)/3;
            Cd(3) = Cd(3)/3;
    end
end
if Cd(end) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
    Cd(end) = (Cd(end-1))/2;
    Cd(end-1)= (Cd(end-1))/2;
    if Cd(end-1) == 0
            Cd(end) = Cd(end-2)/3;
            Cd(end-1) = Cd(end-2)/3;
            Cd(end-2) = Cd(end-2)/3;
    end
end
if Cd(1) < 0
    Cd(1) = -Cd(1);
end
ind0=find(Cd==0);
   if isempty(ind0)==0
    Cd(ind0)=1/3*Cd(ind0-1)+1/3*Cd(ind0+1);
    Cd(ind0-1)=2/3*Cd(ind0-1);
    Cd(ind0+1)=2/3*Cd(ind0+1);
  end


tq = find(QC > 0);
C = Cd(1:end);
[mx,im]=max(Cd); 
tds=td([1:im-2,im+2:length(Cd)]);
ids=[1:im-2,im+2:length(Cd)];
M0=M(tq(1))
Cdd = C;
%Cdd1(1)=Cd1(1);
Td=.1609*C;
Cdd(1)=C(1);
I0g = (C(1)/rho)*T;
% fit values from entire country fit with exception of I0
paramguess = [I0g,(1.66014398889413e-08)*(1e8),1407.12998368299,0.391575667314593];
lb = [I0g*.1,2/T , 0 , .01 ];
ub = [I0g*5,6/T , inf , 1 ];

M = M';
C = C';
 
function X0 = IC(param,C,N)
    %param
    I0 = param(1);
    beta = param(2);
    phi = param(4);
    S0=N;
    
    E0=(beta*(1-(phi*rho))*tau*I0);
    
    Ec0=(beta*(phi*rho)*tau*I0);
    Ic0=phi*rho*I0;
    X0=[S0 E0 I0 (p*tau*M0/(tau+1)) (p*M0/(tau+1)) p*M0 (1-phi)*C(1) phi*rho*C(1) M0];
    %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end


function[z] = Solve_sys(param,td)
    %param
    beta = param(2);
    sigma = param(3);
    phi = param(4);
    X0 = IC(param,C,N);
    %X02(4) = X02(4)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(5) = X02(5)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(6) = X02(6)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(9) = X02(9)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

   [t,Sol] = ode45(@(t,X)Eode(t,X),td,X0);
       
   Cdd(2:end)=Sol(2:end,7)+Sol(2:end,8)-(Sol(1:end-1,7)+Sol(1:end-1,8));
   CumSq=(Sol(tq-3,9))*sqwt;
   CumSq(1) = M(tq(1))*sqwt;
   
   
   Cds = Cdd(ids);
   z=[Cds',CumSq'];
      
      
end


%length(Solve_sys(paramguess,td))
%length([C,M(tq)])
    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C(ids),M(tq)*sqwt],lb,ub);
    [paramfitChina,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C(ids),M(tq)*sqwt],lb,ub)
    
    
    zfit = Solve_sys(paramfitChina,td);
    Cfit = zfit(1:length(tds));
    Qfit = zfit(length(Cfit)+1:end);
    
    
    figure
    plot(tds,Cfit,'b-',td,Cd,'r*','linewidth',2)
    title('China')
    ylabel('New Daily Cases')
    xlabel('Days Since Firt Reoprted Case')
    ylim([0 inf])
    
    figure
    plot(tq,Qfit/sqwt,'b-',tq,M(tq),'r*','linewidth',2)
    title('China')
    ylabel('Quarintined Individuals')
    xlabel('Days Since First Reported Case')
    ylim([0 inf])
    
    %paramfitChina(2) = paramfitChina(2)*T;
   
    %save('paramfitChina.mat','paramfitChina')
%    I0 = paramfitChina(1);
      R0b = paramfitChina(2)*T;
      paramfitChina(2) = R0b;
%    sigma = paramfitChina(3)
%    phi = paramfitChina(4)
    paramfitChina(5) = resnorm;
    
  save('paramfitChina.mat','paramfitChina')
  array2table(paramfitChina)
end

