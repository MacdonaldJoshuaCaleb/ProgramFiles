 function []= Fit_China_simple_model_unreported
filename = '\Users\macdo\Dropbox\COVID-19_modelfitting_data\data-preview_cn.xlsx';
filename = 'data-preview_cn_apr.xlsx';
xlRange = 'E34:BJ37';
dataArray = xlsread(filename,xlRange);
CC=dataArray(1,:)
QC=dataArray(2,:)
M=dataArray(3,:)

% CG=[3,9,8,23,20,15,27,14,14,13,15,8,9,4]; %Beni
% TG=[0,3,4,5,4,3,13,5,3,9,8,3,6,1];
% TGl=[0,2,1,1,1,2,1,3,1,3,4,2,1,0];
% length(CG)
% length(TG)
% length(TGl)
% 
% % theta_v=(TG(2:end)-TGl(2:end))./TG(2:end);
% % theta_v=[1,theta_v];
% % CS=[121,79,84,81,96,63,81,58,55,33,25,9,9,12];
% % TS=[26,32,45,29,18,49,52,39,46,22,14,6,4];
% HG=CG(1:length(TG))-TG;
% % HS=CS(1:length(TS))-TS;
% % p1=exp2cdf(7,10,4);
% % p2=exp2cdf(14,10,4)-p1;
% % p3=1-p1-p2;
Qd=M(2:end);
Cd=diff(CC);
td=1:length(CC);

Qd(isnan(Qd))=0;
QC(isnan(QC))=0;
tq=find(QC>0)
% Q=Qd(Iqp);
% C=Cd(Iqp);
 Q=Qd(1:end);
 C=Cd(1:end);
Nd=length(C);
NQ=length(Q);

C=C';

% [mx,im]=max(Cd)
% csm=sum(Cd(im-1:im+1));
% newavg=mean([Cd(im-2),Cd(im+2:im+5)]);
% Cd(im-1:im+1)=newavg*ones(1,3);
% adcases=csm-3*newavg;
% Cd(1:im)=round(adcases*Cd(1:im)/sum(Cd(1:im-1)))+Cd(1:im);
% Nd=length(Cd);






% 
T=4.64;
T_e=2.71;

 N=1.6563e+08;


tau=3;

sqwt=0.05;
alphac=1/14;
rho=0.5;
 function [dxdt]=Eode(t,X,beta,psi,phi,p,T)
%function [dxdt]=Eode(t,X,beta,psi,phi,p,tau)
        %beta=p(1);
        %T=p(2);
        
       % betaM=thetap*beta;
        %
        %ps=1/sigma*(ddr-(1-p)/p*phi);
        %ps=0;
        %T_e=q(3);
        
     dxdt(1)=-(1+psi)*X(1)*(beta*X(3));
     
     dxdt(2)= beta*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2);
      dxdt(4)= beta*(phi*rho)*X(1)*X(3)- 1/tau*X(4);
     dxdt(3)= 1/tau*X(2)-1/T*X(3);
     dxdt(5)=1/tau*X(4)-1/T_e*X(5);
     dxdt(6)=beta*(phi*rho)*X(1)*X(3);
     dxdt(7)= rho/T*X(3);
     dxdt(8)=1/T_e*X(5);
    % dxdt(9)=beta*(phi)/p*X(1)*X(3)-alphac*X(9);
    dxdt(9)=beta*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5);
     dxdt=dxdt';
 end

 
function[z] = Solve_sys(param,td)
        
       psi=param(3);
       beta=param(4);
      
        p=param(6);
        phi=param(5);
      % tau=param(8);
        
       % T_e=q(3);
        T=param(7);
        I0=param(1);
        
        %Sc0=param(2);
       % S0=param(2);
       Mc0=param(2);
   
%         
% I0=C(1);
% Ic0=Td(1);
%Sc0=Q(1);
        E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N;
        Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
        Ic0=phi*rho*I0;
%         
%         
        S0=N;
%     
     
     X0=[S0 E0 I0 Ec0 Ic0 QC(1) (1-phi)*C(1) phi*C(1) Mc0];
%    
% td0=2:45;
     [t,Sol] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,T),td,X0);
% [t,Sol] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,tau),td,X0);
     
%      figure(6)
%      plot(td0,Sol(:,5)+Sol(:,6),'r',td0,Sol(:,6),'b',1:Nd,C,'r*',1:NQ,Td,'b*')
%       
      CumC=Sol(td,7)+Sol(td,8);
      CumSq=(Sol(tq,4)+Sol(tq,5)+Sol(tq,9))*sqwt;
%       for j=2:length(td)
%           
%               Hh(j-1)=Sol(j,5)-Sol(j-1,5);
%               Tt(j-1)=Sol(j,7)-Sol(j-1,7);
%           
%       end
      z=[CumC',CumSq'];
end
% 
 Tg=4.64;
 %T_eg=2.71;
%phi1=mean(Td(5:14)./C(5:14));

phi1=0.5;
p0=.05;
rho0=.5;
psi0=(1-p0)/p0+600;

betag=2.2807e-08;
%betaMg=c0*p*theta;

 I0g=312/rho0;

 
 Mc0g=phi1*I0g*N*betag;
 paramguess=[I0g,Mc0g,psi0,betag,phi1,p0,Tg];
 lb=zeros(1,length(paramguess));
 lb(4)=betag*.05;
lb(5)=0.05;
 %lb(9)=0.95;
 %lb(4)=T_eg;
 ub=[10*I0g,10*Mc0g,2000,4/N/Tg,1,.2,7];
% 
% td=0:7:Ng*7;
% l1=length(td)
% NG-2
% 
%Q(tq)
 [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[CC,M(tq)*sqwt],lb,ub)
% 
save('parameters_simple_model_unreported_rho50.mat','paramfit');
 psi=paramfit(3)
       beta=paramfit(4)
      %tau=paramfit(8)
        p=paramfit(6)
        phi=paramfit(5)
       
        T=paramfit(7)
       % T_e=q(3);
        
        I0=paramfit(1)
        
        %Sc0=param(2);
        S0=N
Mc0=paramfit(2)
% beta=betaf;
 
 %T_e=Tefit;
 %T=Tfit;
  Rfit=beta*T*S0

% E0=beta*(1-phi)*T*I0/(1/tau*T+1)+betaM*(1-phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1));
       % Ec0=beta*(phi)*T*I0/(1/tau*T+1)+betaM*(phi)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1));
 E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N
        Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
        Ic0=phi*rho*I0;
      
        X0=[S0 E0 I0 Ec0 Ic0 QC(1) (1-phi)*C(1) phi*C(1) Mc0];

        
 Refit=(1-phi*rho)*Rfit

 

% 
 zFit=Solve_sys(paramfit,td);
% 
 CFit=zFit(1:length(td));
 QFit=zFit(length(td)+1:end);
% 
 figure (8)
 plot(td,CFit,'b-',td,CC,'r*')
  figure (9)
 plot(tq,QFit/sqwt,'b-',td,M,'r*')
%  
 [tn,Soln] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,T),td,X0);
%[tn,Soln] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,tau),td,X0);
 
 figure(10)
 plot(tn,Soln(:,7)+Soln(:,8))
 
 St=Soln(:,1);
 It=Soln(:,3);
 Ict=Soln(:,5);

 
 Ret=T*(1-phi)*beta.*St;
 
 
  figure(11)
 plot(tn,Ret)
 





%  ItfigU(1:Nd,[C,Td,Idt,Iut+Idct,Idct])
end