function [] = SimpleModelErlangFitting_C
close all
clear 
warning off;
set(0,'DefaultFigureVisible','on')
% stop/start figures from displaying in matlab, switch 'off' to 'on to display
% import data 
ChinaData = readtable('ChinaDataReshaped.csv');
ProvincePops = readtable('ChinaProvincePops.csv');
ProvincePops = table2array(ProvincePops);
% national total cases
CCt = table2array(ChinaData(:,end-2));
% national total quarintined 
QCt = table2array(ChinaData(:,end));





rho=1;
%sqwt= abs(.01 - .5)/2;
sqwt = .1;
alphac=1/14;
p = .06;
N = sum(ProvincePops); % China pop
ni = 1;
ne = 1;
CC = CCt;

T=4.64;
T_e=2.71;
tau=3;

if ni > 1 || ne > 1
T = T/ni
T_e = T_e/ni
tau = tau/ne
end

td=2:length(CC);
CC = CC(CC > 0);

QC = QCt;
M=table2array(ChinaData(:,end-1));

% provincial daily cases
Cd = diff(CC);
% account for no new reported cases on day after first reported case
if Cd(1) == 0
    Cd(1) = (Cd(2))/2;
    Cd(2)= (Cd(2))/2;
    if Cd(2) == 0
            Cd(1) = Cd(3)/3;
            Cd(2) = Cd(3)/3;
            Cd(3) = Cd(3)/3;
    end
end
if Cd(end) == 0
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
Td=.1609*C;
Cdd(1)=C(1);

    function X0 = IC(param,C,M,tq,N,ni,ne)
        I0 = param(1);
        beta = param(2);
        phi = param(4);
   
        S0=N;
        E0=beta*(1-phi*rho)*tau*I0;
        M0 = M(tq(1));
        temp = [S0 E0 I0 p*tau*M0/(tau+1) p*M0/(tau+1) p*M0 (1-phi)*C(1) phi*rho*C(1) M0];
        
        if ni == 1 && ne == 1
            X0=temp;
        end
        
        if ni > 1 && ne == 1
            X0 = zeros(1,11);
            X0(1) = temp(1);
            X0(2) = temp(2);
            X0(3) = temp(3);
            X0(5) = temp(4);
            X0(6) = temp(5);
            X0(8) = temp(6);
            X0(9)= temp(7);
            X0(10)= temp(8);
            X0(11)= temp(9);
        end
        
        if ni == 1 && ne > 1
            X0 = zeros(1,11);
            X0(1) = temp(1);
            X0(2) = temp(2);
            X0(4) = temp(3);
            X0(5) = temp(4);
            X0(7) = temp(5);
            X0(8) = temp(6);
            X0(9) = temp(7);
            X0(10) = temp(8);
            X0(11) = temp(9);
            
        end
        
        if ni > 1 && ne > 1
            X0 = zeros(1,5+2*ni+2*ne);
            X0(1) = temp(1);
            X0(2) = temp(2);
            X0(4) = temp(3);
            X0(6) = temp(4);
            X0(8) = temp(5);
            X0(10) = temp(6);
            X0(11) = temp(7);
            X0(12) = temp(8);
            X0(13) = temp(9);
        end
    end

 function [dxdt]=EodeErlang(t,X,param,ni,ne,N)
  
    beta = param(2);
    sigma = param(3);
    phi = param(4);
  
    psi = sigma + ((1-p)/p)*phi*rho;

    if ni == 1 && ne == 1
        dxdt = [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';
    end
    
    if ni == 1 && ne > 1
        % S
        dxdt(1) = -(1+psi)*(beta/N)*X(1)*(X(4));
        % Es
        dxdt(2) = (beta/N)*(1-phi*rho)*X(1)*(X(4))-(ne/tau)*X(2);
        dxdt(3) = (ne/tau)*(X(2) - X(3));
        % I
        dxdt(4) = (ne/tau)*X(3)-(ni/T)*X(4); 
        % (E_c)s
        dxdt(5) = (beta/N)*rho*X(1)*(X(4))-(ne/tau)*X(5);
        dxdt(6) = (ne/tau)*(X(5) - X(6));
        % I_c
        dxdt(7) =  (ne/tau)*X(6) - (ni/T_e)*X(7);
        % decoupled terms
        dxdt(8) = (beta/N)*phi*rho*X(1)*(X(4));
        dxdt(9) = (ni/T)*rho*X(4);
        dxdt(10) = (ni/T_e)*rho*X(7);
        dxdt(11) = (beta/N)*((phi*rho)/p)*X(1)*(X(4)) - alphac*X(11);
        
        dxdt = dxdt';
    end
    
    if ni > 1 && ne == 1
        % S =  -(1+psi)*X(1)*((beta/N)*X(3))
        dxdt(1) = -(1+psi)*(beta/N)*X(1)*(X(3) + X(4));
        % E = (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2)
        dxdt(2) = (beta/N)*(1-phi*rho)*X(1)*(X(3) + X(4))-(1/tau)*X(2);
        % Is = 1/tau*X(2)-1/T*X(3)
        dxdt(3) = (1/tau)*X(2)-(ni/T)*X(3);  
        dxdt(4) = (ni/T)*(X(3)-X(4));
       
        % (E_c) = (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4)
        dxdt(5) = (beta/N)*rho*X(1)*(X(3)+X(4))-(1/tau)*X(5);
        % (I_c)s = 1/tau*X(4)-1/T_e*X(5)
        dxdt(6) = (1/tau)*X(5) - (ni/T_e)*X(6);
        dxdt(7) = (ni/T_e)*(X(6) - X(7));
        % decoupled terms
        dxdt(8) = (beta/N)*phi*rho*X(1)*(X(3)+X(4));
        dxdt(9) = (ni/T)*rho*(X(4));
        dxdt(10) = (ni/T_e)*rho*X(7);
        dxdt(11) = (beta/N)*((phi*rho)/p)*X(1)*(X(3)+X(4)) - alphac*X(11);
        
        dxdt = dxdt';
    end
    
    if ni > 1 && ne > 1
        % S
        dxdt(1) = -(1+psi)*(beta/N)*X(1)*(X(4)+X(5));
        % Es
        dxdt(2) = (beta/N)*(1-phi*rho)*X(1)*(X(4)+X(5))-(ne/tau)*X(2);
        dxdt(3) = (ne/tau)*(X(2) - X(3));
        % Is
        dxdt(4) = (ne/tau)*X(3)-(ni/T)*X(4);  
        dxdt(5) = (ni/T)*(X(4) - X(5));
        % (E_c)s
        dxdt(6) = (beta/N)*rho*X(1)*(X(4)+X(5))-(ne/tau)*X(6);
        dxdt(7) = (ne/tau)*(X(6) - X(7));
        % (I_c)s
        dxdt(8) =  (ne/tau)*X(7) - (ne/T_e)*X(8);
        dxdt(9) = (ni/T_e)*(X(8) - X(9));
        % decoupled terms
        dxdt(10) = (beta/N)*phi*rho*X(1)*(X(4)+X(5));
        dxdt(11) = (ni/T)*rho*X(5);
        dxdt(12) = (ni/T_e)*rho*X(9);
        dxdt(13) = (beta/N)*((phi*rho)/p)*X(1)*(X(4)+X(5))- alphac*X(13);
        
        dxdt = dxdt';
    end
 end

 function[z] = Solve_sys(param,td)
    
    X0 = IC(param,C,M,tq,N,ni,ne);

    [t,Sol] = ode45(@(t,X)EodeErlang(t,X,param,ni,ne,N),td,X0);
  
    Cdd(2:end)=Sol(2:end,end-2)+Sol(2:end,end-1)-(Sol(1:end-1,end-2)+Sol(1:end-1,end-1));
    CumSq=(Sol(tq-3,end))*sqwt;
   
   
   CumSq(1) = M(tq(1))*sqwt;
   Cds = Cdd(ids);
   z=[Cds',CumSq'];
      
 end

I0g = (C(1)/rho)*4.64;
% fit values from entire country fit with exception of I0


paramguess = [I0g,(1.66014398889413e-08)*(1e8),1407.12998368299,0.391575667314593];
lb = [I0g*.1,2/(T*ni) , 0 , .01 ];
ub = [I0g*5,6/(T*ni) , inf , 1 ];


 
    
M = M';
C = C';

[paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C(ids),M(tq)*sqwt],lb,ub);
%lb(1) = paramfit(1)*.1;
%ub(1) = paramfit(1)*5;
[paramfitErlang,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C(ids),M(tq)*sqwt],lb,ub);
    
    
    zfit = Solve_sys(paramfitErlang,td);
    Cfit = zfit(1:length(tds));
    Qfit = zfit(length(Cfit)+1:end);
    
    
    figure
    plot(tds,Cfit,'b-',td,Cd,'r*','linewidth',2)
    title('China')
    ylabel('New Daily Cases')
    xlabel('Days Since Firt Reoprted Case')
    
    figure
    plot(tq,Qfit/sqwt,'b-',tq,M(tq),'r*','linewidth',2)
    title('China')
    ylabel('Quarintined Individuals')
    xlabel('Days Since First Reported Case')
    
    
    %paramfitChina(2) = paramfitChina(2)*T;
   
    %save('paramfitChina.mat','paramfitChina')
    I0 = paramfitErlang(1)
   ni*T
      R0b = paramfitErlang(2)*T*ni
      paramfitErlang(2) = R0b;
    sigma = paramfitErlang(3)
    phi = paramfitErlang(4)
    paramfitErlang(5) = resnorm
    
  save('paramfitErlang.mat','paramfitErlang')
array2table(paramfitErlang)
save('tq.mat','tq')
end

