function []= FitSimpleModelAggFixpCI
close all
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
CData2 = table2array(ChinaData(:,15));
CData1 = CCt - CData2;
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
Mc0 = 0.00172837295386061;

%M0=121*exp(r*(Tc(end)-3))/200
N1 = sum(ProvincePops); % China pop
N1 = N1 - ProvincePops(14); % China pop les Hubei 
N2 = ProvincePops(14);
% retrieve province name
%Name = Province(its);
% provincial total cases
CC3 = CData1+CData2;
CC3 = CC3(CC3 > 0);
td3 = 2:length(CC3);
CC1 = CData1;
CC1 = CC1(CC1 > 0)
td1=2:length(CC1);
CC2 = CData2;
CC2 = CC2(CC2 > 0)
%CC2=CC2(2:end);
td2 = 2:length(CC2);

if length(td1) > length(td2)
    td = td1;
end
if length(td1) <= length(td2)
    td = td2;
end
%zers1 = length(CCt)-length(CC1)+1;
%zers2 = length(CCt) - length(CC2)+1;
QC = QCt;
M=table2array(ChinaData(:,end-1));
Qd = M(2:end);
% provincial daily cases
Cd1 = diff(CC1);
% account for no new reported cases on day after first reported case
if Cd1(1) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
    Cd1(1) = Cd1(2);
    if Cd1(1) == 0
            Cd1(1) = Cd1(3);
            Cd1(2) = Cd1(3);
    end
end
if Cd1(1) < 0
    Cd1(1) = -Cd1(1);
end
Cd2 = diff(CC2);
    % account for no new reported cases on day after first reported case
if Cd2(1) == 0
   %Cd(1) = 2*Cd(2)-Cd(3);
    Cd2(1) = Cd2(2);
   % Cd2(2)=0;
    if Cd2(1) == 0
        Cd2(1) = Cd2(3);
        Cd2(2) = Cd2(3);
        Cd2(3)=0;
    end
end
if Cd2(1) < 0
    Cd2(1) = -Cd2(1);
end
Cd3 = diff(CC3);
    % account for no new reported cases on day after first reported case
if Cd3(1) == 0
   %Cd(1) = 2*Cd(2)-Cd(3);
    Cd3(1) = Cd3(2);
   % Cd2(2)=0;
    if Cd3(1) == 0
        Cd3(1) = Cd3(3);
        Cd3(2) = Cd3(3);
        Cd2(3)=0;
    end
end
if Cd3(1) < 0
    Cd3(1) = -Cd3(1);
end

tq = find(QC > 0);
C1 = Cd1(1:end)
C2 = Cd2(1:end)
C3 = Cd3(1:end)
[mx,im]=max(Cd2); 
tds=td2([1:im-2,im+2:length(Cd2)]);
ids=[1:im-2,im+2:length(Cd2)];
[mx,im3]=max(Cd3);
tds3 = td3([1:im3-2,im3+2:length(Cd3)]);
ids3= [1:im3-2,im3+2:length(Cd3)];
C2 = C2([1:im-2,im+2:length(Cd2)]);
C3 = C3([1:im3-2,im3+2:length(Cd3)]);
   M0=M(tq(1))
Cdd1 = C1;
%Cdd1(1)=Cd1(1);
Td1=.1609*C1;
Cdd2 = Cd2;
Cdd2(1)=C2(1);
Td2=.1609*C2;

Cdd3 = Cd3;
Cdd3(1)=C3(1);
Td3=.1609*C3;
    
I01g = (C1(1)/rho)*T;
I02g = (C2(2)/rho)*T;
I03g = (I01g+I02g)*T;
% fit values from entire country fit with exception of I0
paramguess = [I01g,(1.66014398889413e-08)*(1e8),	1407.12998368299,...
              I02g,(1.66014398889413e-08)*(1e8),...
              1407.12998368299,0.391575667314593];
lb = [I01g*.05, 2/T , 0  , I02g*.05, 2/T , 0 , .01 ];
ub = [I01g*10,6/T, inf  ,I02g*10, 6/T , inf , 1 ];

M = M';
C1 = C1';
C2 = C2';
C3 = C3'; 
function X0 = IC(param,QC,C1,C2,C3,M,tq,Td,N,its)
    I0 = param(1);
    beta = param(2);
    if its == 1
        C = C1;
       phi = .217; 
    end
    if its ~= 1
    phi = param(4);
    end
    S0=N;
    
    
    if its == 2
        C = C2;
    
    end
    if its == 3
        C = C3;
    end
    E0=(beta*(1-phi*rho)*tau*I0)*(C(1)/(C1(1)+C2(1)));
    Ec0=(beta*(phi*rho)*tau*I0)*(C(1)/(C1(1)+C2(1)));
    Ic0=phi*rho*I0*(C(1)/(C1(1)+C2(2)));
    X0=[S0 E0 I0 (p*tau*M0/(tau+1))*(C(tq(1))/(C1(tq(1))+C2(tq(1)))) (p*M0/(tau+1))*(C(tq(1))/(C1(tq(1))+C2(tq(1)))) p*M0*(C(tq(1))/(C1(tq(1))+C2(tq(1)))) (1-phi)*C(1) phi*rho*C(1) M0*(C(tq(1))/(C1(tq(1))+C2(tq(1))))];
    %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end




sqwt2 = .01;
function[z] = Solve_sys2(param,td)
    beta = param(2);
    sigma = param(3);
    phi = param(4);
    N = N1+N2;
    param1 = param(1:4);
    X01 = IC(param1,QC,C1,C2,C3,M,tq,Td3,N,1);
    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)+1/T_e*X(5)]';
  

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td1,X01);
   
       
   Cdd3(2:end)=Sol1(2:end,7)+Sol1(2:end,8)-(Sol1(1:end-1,7)+Sol1(1:end-1,8));
   CumSq=(Sol1(tq-3,9))*sqwt2;
   CumSq(1) = M(tq(1))*sqwt2;
   
   
   Cds3 = Cdd3(ids);
   z=[Cds3',CumSq'];
      
      
end

function[z] = Solve_sys(param,td)
    beta = param(2);
    sigma = param(3);
    phi = .217;
    N = N1;
    param1 = param(1:3);
    X01 = IC(param1,QC,C1,C2,C3,M,tq,Td1,N,1);
    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td1,X01);
   
    beta = param(5);
    sigma = param(6);
    phi = param(7);
    N = N2;
    param2 = param(4:7);
    X02 = IC(param2,QC,C1,C2,C3,M,tq,Td2,N,2);
    %X02(4) = X02(4)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(5) = X02(5)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(6) = X02(6)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(9) = X02(9)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

   [t,Sol2] = ode45(@(t,X)Eode(t,X),td2,X02);
       
   Cdd1(2:end)=Sol1(2:end,7)+Sol1(2:end,8)-(Sol1(1:end-1,7)+Sol1(1:end-1,8));
   Cdd2(2:end)=Sol2(2:end,7)+Sol2(2:end,8)-(Sol2(1:end-1,7)+Sol2(1:end-1,8));
   CumSq=(Sol1(tq-3,9))*sqwt + (Sol2(tq-3,9))*sqwt;
   CumSq(1) = M(tq(1))*sqwt;
   Cds1=Cdd1(1:end);
   
   Cds2 = Cdd2(ids);
   z=[Cds1',Cds2',CumSq'];
      
      
end



    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C1,C2,M(tq)*sqwt],lb,ub);
    [paramfitAggFixp,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C1,C2,M(tq)*sqwt],lb,ub)
    I03g = 1500;
    paramguess2 = [I03g,(1.66014398889413e-08)*(1e8),1407.12998368299,0.391575667314593];
    lb2 = [1000,2/T,0,.1];
    ub2 = [2200,6/T,inf,1];

    [paramfit2,resnorm2] = lsqcurvefit(@Solve_sys2,paramguess2,td,[C3,M(tq)*sqwt],lb2,ub2);
    [paramfitAC,resnorm2] = lsqcurvefit(@Solve_sys2,paramfit2,td,[C3,M(tq)*sqwt],lb2,ub2);
    paramfitAC(2) = paramfitAC(2)*T;
    array2table(paramfitAC)
    
    
    % store fit parameters 
    %save('paramfitAgg.mat','paramfitAgg')
    
    zfit = Solve_sys(paramfitAggFixp,td);
    zfit2 = Solve_sys2(paramfitAC,td3);
    Cfit3 = zfit2(1:length(tds3));
    Qfit2 = zfit2(length(Cfit3)+1:end);
    Cfit1 = zfit(1:length(td1));
    Cfit2 = zfit(length(Cfit1)+1:length(tds)+length(Cfit1));
    Qfit = zfit(length(Cfit1)+length(Cfit2)+1:end);
    
    
    figure
    plot(tds3,Cfit3,'b-',td3,Cd3,'r*','linewidth',2)
    title('China')
    ylabel('New Daily Cases')
    xlabel('Days Since Firt Reoprted Case')
    
    figure
    plot(tq,Qfit2/sqwt2,'b-',tq,M(tq),'r*','linewidth',2)
    title('China')
    ylabel('Quarintined Individuals')
    xlabel('Days Since First Reported Case')
    
    figure
    plot(td1,Cfit1,'b-',td1,Cd1,'r*','linewidth',2)
    title('China Less Hubei')
    ylabel('New Daily Cases')
    xlabel('Days Since First Reported Case')
%     baseFileName = sprintf('CasesChinaLessHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    figure
    plot(tds,Cfit2,'b-',td2,Cd2,'r*','linewidth',2)
    title('Hubei')
    ylabel('New Daily Cases')
%     baseFileName = sprintf('CasesHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    figure 
    plot(tq,Qfit/(sqwt),'b-',tq,M(tq),'r*','linewidth',2)
    title('All China')
    xlabel('Days Since First Reported Case')
    ylabel('Quarintined Individuals')
%     baseFileName = sprintf('AggQuarFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    % store fit parameters 
    paramfitAggFixp(2) = paramfitAggFixp(2)*T;
    paramfitAggFixp(5) = paramfitAggFixp(5)*T;
    save('paramfitAggFixp.mat','paramfitAggFixp')
    
      beta = paramfit(2);
    sigma = paramfit(3);
    phi = paramfit(4);
    N = N1;
    param1 = paramfit(1:4);
    X01 = IC(param1,QC,C1,C2,C3,M,tq,Td1,N,1);
    %X01(4) = X01(4)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(5) = X01(5)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(6) = X01(6)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X01(9) = X01(9)*(C1(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td1,X01);
   
   %figure(4)
   %plot(t,Sol1(:,2))
   
    %  figure(5)
   %plot(t,Sol1(:,3))
   
   %figure(6)
   %plot(t,Sol1(:,4))
   
    %  figure(7)
   %plot(t,Sol1(:,5))
   %tq(1)
   
   function [d1]=getData(CFit,nLevel,sets)
    d1 = zeros(sets,length(CFit));
    for k = 1:sets
        for j = 1:length(CFit)
            d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            while d1(k,j) < 0
                d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            end
            d1(k,j) = ceil(d1(k,j));
        end
    end
   end
 %  nLevel1 = .5;
 %  nLevel2 = .5;
 %  sets = 10000;
 %  SynthC1 = getData(Cfit1,nLevel1,sets);
 %  SynthC2 = getData(Cfit2,nLevel2,sets);
   
  % figure
 %  hold on
 %  for w = 1:10000
 %      plot(td1,SynthC1(w,:),'b*-','linewidth',2)
 %  end
 %  plot(td1,Cd1,'r*','linewidth',2)
 %  plot(td1,Cfit1,'k','linewidth',2)
 %  hold off
   
 %  figure
 %  hold on
 %  for w = 1:sets
 %      plot(tds,SynthC2(w,:),'b*-','linewidth',2)
 %  end
 %  plot(td2,Cd2,'r*','linewidth',2)
 %  plot(tds,Cfit2,'k','linewidth',2)
 %  hold off
   
 %  FitparamsSynth = zeros(10000,length(paramfit));
 %  for w = 1:10000
 %      paramguess = paramfitAggFixp;
 %      lb = [paramguess(1)*.05, 2/T , 0, paramguess(4)*.05, 2/T , 0 , .01 ];
 %      ub = [paramguess(1)*10,6/T, inf   ,paramguess(4)*10, 6/T , inf , 1 ];
 %      [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[SynthC1(w,:),SynthC2(w,:),M(tq)*sqwt],lb,ub);
 %      FitparamsSynth(w,:) = paramfit;
 %      if mod(w,100) == 0
 %          fprintf('iteration %i of %i',w,sets)
 %      end
 %  end
 %  save('FitparamsSynth.mat','FitparamsSynth')
end
