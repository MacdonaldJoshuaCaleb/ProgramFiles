function []= FitSimpleModelAggMobility
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

filename = 'WCMI_baidu.xlsx';
xlRange1 = 'U2:AZ363';
xlRange2 = 'U169:AZ184';

dataArray2 = xlsread(filename,xlRange2);

dataArray1 = xlsread(filename,xlRange1);
dataArray1=dataArray1([1:167,184:362],:);
ChAg1=mean(dataArray1);
ChAg2=mean(dataArray2);

[mm1,indd1]=min(ChAg1);
ChAg1=ChAg1(1:indd1);

x=1:length(ChAg1);
%y=Re;
a=ChAg1(1);
    function v=gex(param1,td)
       b=param1(1);
       c=param1(2);
        
       v=(a-b)*exp(-c*(td-1))+b*ones(size(td));
    end


paramg1=[min(ChAg1),-1/length(ChAg1)*log(ChAg1(end)/a)]
% f0=gex(paramg1,x);

[paramfit1,resnorm] = lsqcurvefit(@gex,paramg1,x,ChAg1)

f=gex(paramfit1,x);
figure(1)
plot(x,f,'b',x,ChAg1,'r*')

[mm2,indd2]=min(ChAg2);
ChAg2=ChAg2(1:indd2);

x=1:length(ChAg2);
%y=Re;
a=ChAg2(1);
    


paramg1=[min(ChAg2),-1/length(ChAg2)*log(ChAg2(end)/a)]
% f0=gex(paramg1,x);

[paramfit2,resnorm] = lsqcurvefit(@gex,paramg1,x,ChAg2)

f=gex(paramfit2,x);
figure(2)
plot(x,f,'b',x,ChAg2,'r*')

T=4.64;
T_e=2.71;
tau=3;
rho=1;
sqwt=0.01;
alphac=1/14;
p = .06;
Mc0 = 0.00172837295386061;
N1 = sum(ProvincePops); % China pop
N1 = N1 - ProvincePops(14); % China pop les Hubei 
N2 = ProvincePops(14);
% retrieve province name
%Name = Province(its);
% provincial total cases
CC1 = CData1;
CC1 = CC1(CC1 > 0);
td1=2:length(CC1);
CC2 = CData2;
CC2 = CC2(CC2 > 0);
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
    Cd1(1) = (Cd1(2))/2;
    Cd1(2)= (Cd1(2))/2;
    if Cd1(2) == 0
            Cd1(1) = Cd1(3)/3;
            Cd1(2) = Cd1(3)/3;
            Cd1(3) = Cd1(3)/3;
    end
end
if Cd1(end) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
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
    % account for no new reported cases on day after first reported case
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
paramguess = [I01g,(1.66014398889413e-08)*(1e8),	1407.12998368299,...
              0.391575667314593,I02g,(1.66014398889413e-08)*(1e8),...
              1407.12998368299,0.391575667314593];
lb = [I01g*.5, 2/T , 0 , .01 , I02g*.5, 2/T , 0 , .01];
ub = [I01g*2,6/T, inf , 1  ,I02g*2, 6/T , inf , 1];

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
    %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end

function X0 = IC2(param,QC,C1,C2,M,tq,Td,N,its)
    I0 = param(1);
    beta = param(2);
    phi = param(4);
    p = param(5);
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
    %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end



function[z] = Solve_sys(param,td)
    beta = param(2);
    sigma = param(3);
    phi = param(4);
    N = N1;
    param1 = param(1:4);
    X01 = IC(param1,QC,C1,C2,M,tq,Td1,N,1);
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
   
    beta = param(6);
    sigma = param(7);
    phi = param(8);
    N = N2;
    param2 = param(5:8);
    X02 = IC(param2,QC,C1,C2,M,tq,Td2,N,2);
    %X02(4) = X02(4)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(5) = X02(5)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(6) = X02(6)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(9) = X02(9)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5)]';

   [t,Sol2] = ode45(@(t,X)Eode(t,X),td2,X02);
   Cdd1(2:end)=Sol1(2:end,7)+Sol1(2:end,8)-(Sol1(1:end-1,7)+Sol1(1:end-1,8));
   Cdd2(2:end)=Sol2(2:end,7)+Sol2(2:end,8)-(Sol2(1:end-1,7)+Sol2(1:end-1,8));
   CumSq=(Sol1(tq-3,9))*sqwt + (Sol2(tq-3,9))*sqwt;
  % CumSq(1) = M(tq(1))*sqwt;
   Cds1=Cdd1(1:end);
   Cds2 = Cdd2(ids);
   z=[Cds1',Cds2',CumSq'];
      
      
end
    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C1,C2,M(tq)*sqwt],lb,ub);
    [paramfitAgg,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C1,C2,M(tq)*sqwt],lb,ub);
    % store fit parameters 
    %save('paramfitAgg.mat','paramfitAgg')
     beta = paramfitAgg(2);
    sigma = paramfitAgg(3);
    phi = paramfitAgg(4);
    N = N1;
    param1 = paramfitAgg(1:4);
    X01 = IC(param1,QC,C1,C2,M,tq,Td1,N,1);
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
   
   figure
   plotyy(t,Sol1(:,6)/phi*psi,1:length(ChAg1),ChAg1)
   
    beta = paramfitAgg(6);
    sigma = paramfitAgg(7);
    phi = paramfitAgg(8);
    N = N2;
    param2 = paramfitAgg(5:8);
    X02 = IC(param2,QC,C1,C2,M,tq,Td2,N,2);
    %X02(4) = X02(4)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(5) = X02(5)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(6) = X02(6)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    %X02(9) = X02(9)*(C2(tq(1))/(C1(tq(1))+C2(tq(1))));
    psi = sigma + ((1-p)/p)*phi*rho;
    Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5)]';

   [t2,Sol2] = ode45(@(t,X)Eode(t,X),td2,X02);
   
 figure
   plotyy(t2,Sol2(:,6)/phi*psi,1:length(ChAg2),ChAg2)
   
   
    
    zfit = Solve_sys(paramfitAgg,td);
    Cfit1 = zfit(1:length(td1));
    Cfit2 = zfit(length(Cfit1)+1:length(tds)+length(Cfit1));
    Qfit = zfit(length(Cfit1)+length(Cfit2)+1:end);
    
    
    
 %   figure
 %   plot(td1,Cfit1,'b-',td1,Cd1,'r*','linewidth',2)
 %   title('China Less Hubei')
 %   ylabel('New Daily Cases')
 %   xlabel('Days Since First Reported Case')
%     baseFileName = sprintf('CasesChinaLessHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
%    figure
%    plot(tds,Cfit2,'b-',td2,Cd2,'r*','linewidth',2)
%    title('Hubei')
%    ylabel('New Daily Cases')
%     baseFileName = sprintf('CasesHubeiFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
%    figure 
%    plot(tq,Qfit/(sqwt),'b-',tq,M(tq),'r*','linewidth',2)
%    title('All China')
%    xlabel('Days Since First Reported Case')
 %   ylabel('Quarintined Individuals')
%     baseFileName = sprintf('AggQuarFixp.png');
%     fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\AggPlots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
    % store fit parameters 
%     paramfitAggFixp(2) = paramfitAggFixp(2)*T;
%     paramfitAggFixp(6) = paramfitAggFixp(6)*T;
%     paramfitAggFixp(9) = resnorm;
%     save('paramfitAggFixp_c2.mat','paramfitAggFixp')
%     array2table(paramfitAggFixp)

function dydt=Edde1(t,X,beta,psic,N,phi,p)
        if t<=length(ChAg1)
            sigmat=paramfit1(2);
        else sigmat=0;
        end
        dydt=[-(1+psic)*X(1)*((beta/N)*X(3))-sigmat*X(1),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';
end 

function dydt=Edde2(t,X,beta,psic,N,phi,p)
        if t<=length(ChAg2)
            sigmat=paramfit2(2);
        else sigmat=0;
        end
        dydt=[-(1+psic)*X(1)*((beta/N)*X(3))-sigmat*X(1),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';
end 

clear p;

function [z] = Solve_sys2(param,td)
    beta = param(2);
  % sigma = paramfitAgg(3);
    phi = param(3);
    p = param(4);
    N = N1;
    psic =  ((1-p)/p)*phi*rho;
    param1 = [param(1:2),paramfit1(2),param(3),param(4)];
    X01 = IC2(param1,QC,C1,C2,M,tq,Td1,N,1);

   [ts,Sol1] = ode45(@(t,X)Edde1(t,X,beta,psic,N,phi,p),td,X01);
       
   Cdd1(2:length(td))=Sol1(2:end,7)+Sol1(2:end,8)-(Sol1(1:end-1,7)+Sol1(1:end-1,8));
    
   
    
    beta = param(6);
    phi = param(7);
    p = param(8);
    psic =  ((1-p)/p)*phi*rho;
    N = N2;
    param2 = [param(5:6),paramfit2(2),param(7),param(8)];
    X02 = IC2(param2,QC,C1,C2,M,tq,Td2,N,2);
    
    [ts,Sol] = ode45(@(t,X)Edde2(t,X,beta,psic,N,phi,p),td,X02);
       
   Cdd2(2:length(td))=Sol2(2:end,7)+Sol2(2:end,8)-(Sol2(1:end-1,7)+Sol2(1:end-1,8));
   CumSq=(Sol1(tq-3,9))*sqwt + (Sol2(tq-3,9))*sqwt;
    Cds1=Cdd1(1:end);
   Cds2 = Cdd2(ids);
   z=[Cds1',Cds2',CumSq'];
end
 
temp = paramfitAgg;
temp(3) = [];
temp(6) = [];
paramguess = [temp(1:3),.06,temp(4:6),.06];

templb = lb;
tempub = ub;
templb(3) = [];
templb(6) = [];
tempub(3) = [];
tempub(6) = [];
lb = [templb(1:3),.01,templb(4:6),.01];
ub = [tempub(1:3),.1,tempub(4:6),.1];

[paramfit,resnorm] = lsqcurvefit(@Solve_sys2,paramguess,td,[C1,C2,M(tq)*sqwt],lb,ub);
lb(1) = .5*paramfit(1);
lb(5) = .5*paramfit(5);
ub(1) = 2*paramfit(1);
ub(5) = 2*paramfit(5);
[paramfit,resnorm] = lsqcurvefit(@Solve_sys2,paramfit,td,[C1,C2,M(tq)*sqwt],lb,ub);
 zfit = Solve_sys2(paramfit,td);
    Cfit1 = zfit(1:length(td1));
    Cfit2 = zfit(length(Cfit1)+1:length(tds)+length(Cfit1));
    Qfit = zfit(length(Cfit1)+length(Cfit2)+1:end);
 
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
 paramfit(2) = paramfit(2)*T;
 paramfit(6) = paramfit(6)*T;
array2table(paramfit)
paramfitMobility = paramfit;
save('paramfitMobility.mat','paramfitMobility')
end
