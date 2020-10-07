function []= FitSimpleModelSerialIntVar

%close all
warning off;
set(0,'DefaultFigureVisible','on')
% stop/start figures from displaying in matlab, switch 'off' to 'on to display
%set(0,'DefaultFigureVisible','off')
% import data 
ChinaData = readtable('ChinaDataReshaped.csv');
ProvincePops = readtable('ChinaProvincePops.csv');
ProvincePops = table2array(ProvincePops);
% national total cases
CCt = table2array(ChinaData(:,end-2));
CData2 = table2array(ChinaData(:,15));
CData1 = CCt - CData2;
% national total quarintined 
QCt = table2array(ChinaData(:,end));


T=4.64;
%T=6;
T_e=2.71;
tau=3;
rho=1;
sqwt=0.01;
alphac=1/14;
p = .06;

N1 = sum(ProvincePops); % China pop
N1 = N1 - ProvincePops(14); % China pop les Hubei 
N2 = ProvincePops(14);

CC1 = CData1;
CC1 = CC1(CC1 > 0);
td1=2:length(CC1);
CC2 = CData2;
CC2 = CC2(CC2 > 0)

td2 = 2:length(CC2);
if length(td1) > length(td2)
    td = td1;
end
if length(td1) <= length(td2)
    td = td2;
end

QC = QCt;
M=table2array(ChinaData(:,end-1));
Qd = M(2:end);
% provincial daily cases
Cd1 = diff(CC1);
% account for data irregularities 
if Cd1(1) == 0
    Cd1(1) = (Cd1(2))/2;
    Cd1(2)= (Cd1(2))/2;
    if Cd1(2) == 0
            Cd1(1) = Cd1(3)/3;
            Cd1(2) = Cd1(3)/3;
            Cd1(3) = Cd1(3)/3;
    end
end
if Cd1(end) == 0
    Cd1(end) = (Cd1(end-1))/2;
    Cd1(end-1)= (Cd1(end-1))/2;
    if Cd1(end-1) == 0
            Cd1(end) = Cd1(end-2)/3;
            Cd1(end-1) = Cd1(end-2)/3;
            Cd1(end-2) = Cd1(end-2)/3;
    end
end
if Cd1(1) < 0
    Cd1(1) = -Cd1(1);
end
ind0=find(Cd1==0);
   if isempty(ind0)==0
    Cd1(ind0)=1/3*Cd1(ind0-1)+1/3*Cd1(ind0+1);
    Cd1(ind0-1)=2/3*Cd1(ind0-1);
    Cd1(ind0+1)=2/3*Cd1(ind0+1);
   end
    
Cd2 = diff(CC2);
if Cd2(1) == 0
    Cd2(1) = (Cd2(2))/2;
    Cd2(2)= (Cd2(2))/2;
    if Cd2(2) == 0
            Cd2(1) = Cd2(3)/3;
            Cd2(2) = Cd2(3)/3;
            Cd2(3) = Cd2(3)/3;
    end
end
if Cd2(end) == 0
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
tq = find(QC > 0);
C1 = Cd1(1:end);
C2 = Cd2(1:end);
[mx,im]=max(Cd2); 
ids=[1:im-2,im+2:length(Cd2)];
tds=td(ids);
C2 = C2(ids);
   
Cdd1 = C1;
Cdd1(1)=Cd1(1);
Td1=.1609*C1;
Cdd2 = Cd2;
Cdd2(1)=Cd2(1);
Td2=.1609*C2;

I01g = (C1(1)/rho)*T;
I02g = (C2(2)/rho)*T;
% fit values from entire country fit with exception of I0
paramguess = [632.864026481429,5.99999999986307,74367.0214604314,...
              0.308064520308229,1018.50380516702,5.99999999999980,...
              561.569333839141,0.0489799044337563,9923332.46179674]
lb = [632.864026481429*.5, 2/T , 79769*.5 , .01 , 1018.50380516702*.25, 2/T , 0 , .01];
ub = [632.864026481429*2,6/T, 79769*2 , 1  ,1018.50380516702*2, 6/T , inf , 1];

M = M';
C1 = C1';
C2 = C2';
    
function X0 = IC(param,QC,C1,C2,M,tq,Td,N,its)
    I0 = param(1);
    beta = param(2);
    phi = param(4);
    S0=N;
    
    if its == 1
        C = C1;
        
    end
    if its == 2
        C = C2;
    
    end
    E0=beta*(1-phi*rho)*tau*I0;
    Ec0=beta*(phi*rho)*tau*I0;
    Ic0=phi*rho*I0;
    X0=[S0 E0 I0 Ec0 Ic0 p*QC(1) C(1) phi*rho*C(1) M(tq(1))*(C(tq(1))/(C1(tq(1))+C2(tq(1))))];        
end



function[z] = Solve_sys(param,td)
    beta = param(2);
    sigma = param(3);
    phi = param(4);
    N = N1;
    param1 = param(1:4);
    X01 = IC(param1,QC,C1,C2,M,tq,Td1,N,1);

    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),...
                   (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2) - (sigma*beta)*(X(3)/N)*X(2),...
                   1/tau*X(2)-1/T*X(3) - (sigma*beta)*(X(3)/N)*X(3),...
                   (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                   1/tau*X(4)-1/T_e*X(5),...
                   (beta/N)*(phi*rho)*X(1)*X(3),...
                   rho/T*X(3),...
                   1/T_e*X(5),...
                  (beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td1,X01);
   
    beta = param(6);
    sigma = param(7);
    phi = param(8);
    N = N2;
    param2 = param(5:8);
    X02 = IC(param2,QC,C1,C2,M,tq,Td2,N,2);

    psi = sigma + ((1-p)/p)*phi*rho;
     Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),...
                   (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2) - (sigma*beta)*(X(3)/N)*X(2),...
                   1/tau*X(2)-1/T*X(3) - (sigma*beta)*(X(3)/N)*X(3),...
                   (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                   1/tau*X(4)-1/T_e*X(5),...
                   (beta/N)*(phi*rho)*X(1)*X(3),...
                   rho/T*X(3),...
                   1/T_e*X(5),...
                  (beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

              
   [t,Sol2] = ode45(@(t,X)Eode(t,X),td2,X02);
       
   Cdd1(2:end)=Sol1(2:end,7)+Sol1(2:end,8)-(Sol1(1:end-1,7)+Sol1(1:end-1,8));
   Cdd2(2:end)=Sol2(2:end,7)+Sol2(2:end,8)-(Sol2(1:end-1,7)+Sol2(1:end-1,8));
   CumSq=(Sol1(tq-(tq(1)-1),9))*sqwt + (Sol2(tq-3,9))*sqwt;
   

   Cds1=Cdd1(1:end);
   Cds2 = Cdd2(ids);
   z=[Cds1',Cds2',CumSq'];
      
      
end
    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C1,C2,M(tq)*sqwt],lb,ub);
    
    [paramfitAggFixp,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C1,C2,M(tq)*sqwt],lb,ub);
    % store fit parameters 
    %save('paramfitAgg.mat','paramfitAgg')
    
    zfit = Solve_sys(paramfitAggFixp,td);
    Cfit1 = zfit(1:length(td1));
    Cfit2 = zfit(length(Cfit1)+1:length(tds)+length(Cfit1));
    Qfit = zfit(length(Cfit1)+length(Cfit2)+1:end);
    
    figure(1)
    plot(td1,Cfit1,'b-',td1,Cd1,'r*','linewidth',2)
    title('China Less Hubei')
    ylabel('New Daily Cases')
    xlabel('Days Since First Reported Case')
%     baseFileName = sprintf('CasesChinaLessHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    figure(2)
    plot(tds,Cfit2,'b-',td2,Cd2,'r*','linewidth',2)
    title('Hubei')
    ylabel('New Daily Cases')
    xlabel('Days Since First Reported Case')
%     baseFileName = sprintf('CasesHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    figure(3)
    plot(tq,Qfit/(sqwt),'b-',tq,M(tq),'r*','linewidth',2)
    title('All China')
    xlabel('Days Since First Reported Case')
    ylabel('Quarintined Individuals')
%     baseFileName = sprintf('AggQuarFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    % store fit parameters 

 beta = paramfitAggFixp(2);
    sigma = paramfitAggFixp(3);
    phi = paramfitAggFixp(4);
     Re=beta*T*(1-phi);
    N = N1;
    param1 = paramfitAggFixp(1:4);
    X01 = IC(param1,QC,C1,C2,M,tq,Td1,N,1);
    X01(10)=0;
    psi = sigma + ((1-p)/p)*phi*rho;
Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),...
                   (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2) - (sigma*beta)*(X(3)/N)*X(2),...
                   1/tau*X(2)-1/T*X(3) - (sigma*beta)*(X(3)/N)*X(3),...
                   (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                   1/tau*X(4)-1/T_e*X(5),...
                   (beta/N)*(phi*rho)*X(1)*X(3),...
                   rho/T*X(3),...
                   1/T_e*X(5),...
                  (beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*beta/N*X(3)*(X(1)+X(2)+X(3))]';

   [t,Sol1n] = ode45(@(t,X)Eode(t,X),td1,X01);
   
    beta = paramfitAggFixp(6);
    sigma = paramfitAggFixp(7);
    phi = paramfitAggFixp(8);
     Re=beta*T*(1-phi);
    N = N2;
    param2 = paramfitAggFixp(5:8);
    X02 = IC(param2,QC,C1,C2,M,tq,Td2,N,2);
      X02(10)=0;
    psi = sigma + ((1-p)/p)*phi*rho;
Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),...
                   (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2) - (sigma*beta)*(X(3)/N)*X(2),...
                   1/tau*X(2)-1/T*X(3) - (sigma*beta)*(X(3)/N)*X(3),...
                   (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                   1/tau*X(4)-1/T_e*X(5),...
                   (beta/N)*(phi*rho)*X(1)*X(3),...
                   rho/T*X(3),...
                   1/T_e*X(5),...
                  (beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*beta/N*X(3)*(X(1)+X(2)+X(3))]';

   [t2,Sol2n] = ode45(@(t,X)Eode(t,X),td2,X02);
     
    paramfitAggFixp(2) = paramfitAggFixp(2)*T;
    paramfitAggFixp(6) = paramfitAggFixp(6)*T;
    paramfitAggFixp(9) = resnorm;
   
    
    
    save('paramfitSerialIntVar.mat','paramfitAggFixp')
   
    load('paramfitAgg.mat','paramfitAgg')
     beta = paramfitAgg(2)/T;
    sigma = paramfitAgg(3);
    phi = paramfitAgg(4);
     Re=beta*T*(1-phi);
    N = N1;
    param1 = paramfitAgg(1:4);
    X01 = IC(param1,QC,C1,C2,M,tq,Td1,N,1);
      X01(10)=0;
    psi = sigma + ((1-p)/p)*phi*rho;
     Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*beta/N*X(3)*(X(1))]';

   [t,Sol1] = ode45(@(t,X)Eode(t,X),td1,X01);
   
    beta = paramfitAgg(6)/T;
    sigma = paramfitAgg(7);
    phi = paramfitAgg(8);
     Re=beta*T*(1-phi);
    N = N2;
    param2 = paramfitAgg(5:8);
    X02 = IC(param2,QC,C1,C2,M,tq,Td2,N,2);
      X02(10)=0;
    psi = sigma + ((1-p)/p)*phi*rho;
      Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9),sigma*beta/N*X(3)*(X(1))]';

   [t2,Sol2] = ode45(@(t,X)Eode(t,X),td2,X02);
   
   figure(4)
   plot(t,Sol1n(:,10),'r-',t,Sol1(:,10),'b-')
     
   figure(5)
   plot(t2,Sol2n(:,10),'r-',t2,Sol2(:,10),'b-')
    
   array2table(paramfitAggFixp)
end