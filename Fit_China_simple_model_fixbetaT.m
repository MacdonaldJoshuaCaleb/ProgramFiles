 function []= Fit_avg_Re_exp
% 
% filename = '/Users/cxb0559/Dropbox/Job_applications/nCoV_model_contact_tracing/MATLAB_CODES/data-preview_cn.xlsx';
filename = 'data-preview_cn_apr.xlsx';
xlRange = 'E34:BJ37';
dataArray = xlsread(filename,xlRange);
CC=dataArray(1,:)
QC=dataArray(2,:)
M=dataArray(3,:)
R=dataArray(4,:)
ddr=QC(end)/CC(end);
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
% [mx,im]=max(Cd)
% csm=sum(Cd(im-1:im+1));
% newavg=mean([Cd(im-2),Cd(im+2:im+5)]);
% Cd(im-1:im+1)=newavg*ones(1,3);
% adcases=csm-3*newavg;
% Cd(1:im)=round(adcases*Cd(1:im)/sum(Cd(1:im-1)))+Cd(1:im);
% Nd=length(Cd);

siran=(1:20);
Ns=length(siran);
ctran=(1:12);
% wbscv=repmat([6.72],size(siran));
% wbshv=repmat([5.81],size(siran));
%sicd=wblcdf(siran,wbscv,wbshv);

%Iqp=find(Qd>0);
Qd(isnan(Qd))=0;
QC(isnan(QC))=0;
tq=find(QC>0)
% Q=Qd(Iqp);
% C=Cd(Iqp);
 Q=Qd(1:end);
 C=Cd(1:end);
Nd=length(C);
NQ=length(Q);
NT=NQ-length(ctran);
Nr=NQ-length(siran);
Nc=Nd-length(siran);
Q=Q';
C=C';
%Cd=Cd';
Td=zeros(size(Q));
Re=zeros(Nr,1);
R0=zeros(Nc,1);
Reno=zeros(Nr,1);
Renop=Reno;

    function u=CTicdf(tv)
%         tv=flipud(tv);
%       u=0.75*conv(lognpdf(tv,1.57,0.65),logncdf(tv,0.74,0.64));
%        u=u(1:length(tv))+0.25*logncdf(tv,1.57,0.65);
 u=0.75*logncdf(tv,1.57,0.65);
      u(tv==0)=u(tv==0)+.25*u(tv==0);
     u=flipud(u);   
    end
    function u=sicdCT(tv)
      u=gamcdf(tv,1.8,1/.5);
        
    end
    function u=sicd(tv)
      u=gamcdf(tv,2.29,1/.36);
        
    end
Ncc=length(CC);
figure(7)
plot(1:Ncc,CC)

p=.05;  %prob. of infection/contact
theta=0.16;
% rho=0.5;


for i=1:length(ctran)
    Td(i)=p*(CTicdf(1:i)-CTicdf((1:i)-1))*Q(1:i); 
    while Td(i)>=C(i)
      Td(i)=0.5*Td(i);
     % jj=jj+1;
    end  %traced infected from quarantined and traced to infected cdf%traced infected from quarantined and traced to infected cdf
end

for i=length(ctran)+1:NQ
%     CTicdf(ctran)-CTicdf(ctran-1)
    Td(i)=p*(CTicdf(ctran)-CTicdf(ctran-1))*Q(i-ctran); 
   % jj=1;
%    if Td(i)>=C(i)
%      %uu=Td(i);
%       Td(i)=C(i)*mean(Td(1:i-1)./C(1:i-1));
%       
%      % jj=jj+1;
%     end
    while Td(i)>=C(i)
        uu=Td(i);
      Td(i)=0.5*uu;
     % jj=jj+1;
    end%traced infected from quarantined and traced to infected cdf
end

Td=.1609*C;
H=C-Td;

for i=1:Nr
    Re(i)=1/(H(i)+theta*Td(i))*((sicd(siran)-sicd(siran-1))*H(i+siran)+theta*(sicdCT(siran)-sicdCT(siran-1))*Td(i+siran)); 
    Reno(i)=1/(H(i)+theta*Td(i))*((sicd(siran)-sicd(siran-1))*H(i+siran)+(sicdCT(siran)-sicdCT(siran-1))*Td(i+siran));%traced infected from quarantined and traced to infected cdf
   % Renop(i)=1/(H(i))*((sicd(siran)-sicd(siran-1))*(H(i+siran).*(H(i+siran)+rho*Td(i+siran))./(H(i+siran)+(1-rho)*Td(i+siran))));
end

for i=1:Nc
    R0(i)=1/(C(i))*((sicd(siran)-sicd(siran-1))*C(i+siran)); 
end

% figure(1)
% plot(1:Nr,Re,'b',1:Nr,Reno,'r',1:Nc,R0,'k')
% 
% figure(2)
% plot(1:Nd,C,'r-*',1:NQ,Td,'b-*')

mr1=mean(Re(1:9));
mr2=mean(Re(10:end));


figure(3)
plot(1:Nr,Re,'b',1:Nr,Reno,'r')
%hold on
%plot([1,9],[mr1,mr1],'k--',[10,Nr],[mr2,mr2],'k--')

x=(1:Nr)';
y=Re;
a=Re(1);
    function v=gex(param1,td)
       b=param1(1);
       c=param1(2);
        
       v=(a-b)*exp(-c*td)+b*ones(size(td));
    end
size(Re)

paramg1=[min(Re),-1/Nr*log(Re(end)/a)]
f0=gex(paramg1,x)
size(f0)
[paramfit,resnorm] = lsqcurvefit(@gex,paramg1,x,Re)

f=gex(paramfit,x);

figure(4)
plot(x,f,'b',x,Re,'r*')

figure(5)
plot(1:Nd,Td./C)
% 
% p1=wblcdf(7,6.72,5.81);
% p2=wblcdf(14,13.2,2.6)-p1;
% p3=1-p1-p2;

% NG=length(TG)-1;
% % NS=length(TS)-1;
% R0G=zeros(NG-2,1);
% % R0S=zeros(NS-2,1);
% R0Gno=zeros(NG-2,1);
% % R0Sno=zeros(NS-2,1);
% 
T=4.64;
T_e=2.71;
Tg=T;
T_eg=T_e;
 N=1e+08;
 betaM=0;
% 
% Ng=NG+1;
% % Ns=NS+1;
%  theta=0.13;
% sigma=1/10;
% for i=1:NG
% R0G(i)=(p3*(CG(i+3)-(1-theta_v(i+3))*TG(i+3))+p2*(CG(i+2)-(1-theta_v(i+2))*TG(i+2))+p1*(CG(i+1)-(1-theta_v(i+1))*TG(i+1)))/(CG(i)-(1-theta_v(i))*TG(i));
% R0Gno(i)=(p3*(CG(i+3))+p2*(CG(i+2))+p1*(CG(i+1)))/(CG(i)-(1-theta_v(i))*TG(i));
% end
% % 
% % for i=1:NS-2
% % R0S(i)=(p3*(CS(i+3)-(1-theta)*TS(i+3))+p2*(CS(i+2)-(1-theta)*TS(i+2))+p1*(CS(i+1)-(1-theta)*TS(i+1)))/(CS(i)-(1-theta)*TS(i));
% % R0Sno(i)=(p3*(CS(i+3))+p2*(CS(i+2))+p1*(CS(i+1)))/(CS(i)-(1-theta)*TS(i));
% % end
% 
% qG=TG./CG(1:length(TG));
% 
% qGc=(TG-TGl)./CG(1:length(TG));
% % qS=TS./CS(1:length(TS));
% 
% GavgR=mean(R0Gno)
% %SavgR=mean(R0Sno)
% GavgRe=mean(R0G)
% avgqG=mean(qG)
% %avgqS=mean(qS)
% % 
tau=3;
rho=1;
%tau=3.95;
% Tg=4.64;
% T_eg=2.71;
% %phi1=mean(Td(5:14)./C(5:14));
% phi2=mean(Td(15:end)./C(15:end));
% phi1=0.5;
% sigma0=(1-p)/p
% c0=mean(Re(1:4))/(p*(1-phi1)*Tg+p*phi1*theta*Tg/T_eg*T_eg)/N;
% betag=c0*p;
% betaMg=c0*p*theta;
% Re1=(betag*(1-phi1)*Tg+betaMg*phi1*T_eg)*N
alpha0=1/14;
nu0=0;
 Sc0g=700;
sqwt=0.05;
alphac=1/14;

beta=6/N/T;
 function [dxdt]=Eode(t,X,beta,psi,phi,p,rho,T)
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

 function [dxdt]=EodeL(t,X,beta,phi,N)
%function [dxdt]=Eode(t,X,beta,psi,phi,p,tau)
        %beta=p(1);
        %T=p(2);
        
       % betaM=thetap*beta;
        %
        %ps=1/sigma*(ddr-(1-p)/p*phi);
        %ps=0;
        %T_e=q(3);
        
     
     
     dxdt(1)= beta*(1-phi)*N*X(2)- 1/tau*X(1);
      
     dxdt(2)= 1/tau*X(1)-1/T*X(2);
     
     dxdt=dxdt';
 end
function[y] = Solve_sysL(par,tl)
    E00=par(1);
    I00=par(2);
    [tlb,Sollb] = ode45(@(t,X)EodeL(t,X,beta,phi,N,rho),tl,[E00,I00]);
    y=Sollb(end,:);
end
% 
function[z] = Solve_sys(param,td)
        
       psi=param(3);
     %  beta=param(4);
      
        p=param(5);
        phi=param(4);
      % tau=param(8);
        
       % T_e=q(3);
       % rho=param(6);
        I0=param(1);
      %  T=param(8);
        %Sc0=param(2);
       % S0=param(6);
       Mc0=param(2);
   S0=N;
 %  beta=3/N/T;
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
     X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
    % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) phi*C(1) Mc0];
%    
% td0=2:45;
     [t,Sol] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,rho,T),td,X0);
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
% Tg=4.64;
 %T_eg=2.71;
%phi1=mean(Td(5:14)./C(5:14));
phi2=mean(Td(15:end)./C(15:end));
phi1=0.5;
p0=.05;
rho0=.5;
psi0=(1-p0)/p0+600;
c0=mean(Re(1:4))/((1-phi1)*Tg)/N;
betag=beta;
%betaMg=c0*p*theta;
Re1=(betag*(1-phi1)*Tg)*N
 I0g=312/rho0;

 
 Mc0g=phi1*I0g*N*betag;
 paramguess=[I0g,Mc0g,psi0,phi1,p0];
 lb=zeros(1,length(paramguess));
% lb(4)=betag*.05;
lb(4)=0.05;
 %lb(9)=0.95;
 %lb(4)=T_eg;
 ub=[10*I0g,10*Mc0g,2000,1,.1];
% 
% td=0:7:Ng*7;
% l1=length(td)
% NG-2
% 
%Q(tq)
 [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[CC,M(tq)*sqwt],lb,ub)
% 
save('parameters_simple_lesst.mat','paramfit');
 psi=paramfit(3)
      % beta=paramfit(4)
      %T=paramfit(8)
        p=paramfit(5)
        phi=paramfit(4)
       
     %   rho=paramfit(6)
       % T_e=q(3);
        
        I0=paramfit(1)
        
        %Sc0=param(2);
        %S0=paramfit(6)
        S0=N;
      %  beta=4/N/T
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
        X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
       % X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) phi*C(1) Mc0];
%  A=[-1/tau,beta*S0*(1-phi);1/tau,-1/T]
%  FS=expm(A)
%  FSI=inv(FS)
%  x00=(FSI^40)*[E0;I0]
%  
%  E00=x00(1)
%  I00=x00(2)
 %Rmfit=betaf*T_e*thetafit*S0
 Refit=(1-phi*rho)*Rfit
% betaf-betag
% phif-phig
 fracCT=1
% betaf-betag
% phif-phig
% 
 zFit=Solve_sys(paramfit,td);
% 
 CFit=zFit(1:length(td));
 QFit=zFit(length(td)+1:end);
% 
% cumquarfig(tq,QFit/sqwt,td,M)
 figure (8)
  plot(td,CFit,'b-',td,CC,'r*')
  figure (9)
 plot(tq,QFit/sqwt,'b-',td,M,'r*')
%  
 [tn,Soln] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,rho,T),td,X0);
%[tn,Soln] = ode45(@(t,X)Eode(t,X,beta,psi,phi,p,tau),td,X0);
 
 figure(10)
 plot(tn,Soln(:,7)+Soln(:,8))
 
 St=Soln(:,1);
 It=Soln(:,3);
 Ict=Soln(:,5);
 Itt=It+Ict;
 
 pm=Ict./ (It+Ict);
 pmc=diff(Soln(:,8))./(diff(Soln(:,7))+diff(Soln(:,8)));
 Idt=(diff(Soln(:,7))+diff(Soln(:,8)));
% Iut=(diff(Soln(:,3))+diff(Soln(:,5)));
  Iut=diff(Soln(:,7))/rho;
 Idct=diff(Soln(:,8));
 Lam=2*(1-pm.^2).*(1-phi)*beta.*St.*It./(Itt.^2);
 eps=2*Lam(1)*94;
 Ne=eps./Lam;
 Ret=T*(1-phi)*beta.*St;
 figure(13)
 plot(tn,Ne)
 
  figure(11)
 plot(tn,Ret)
 
%  paramguessL=[.2,.1];
%  lb=zeros(1,length(paramguess));
% 
% 
%  %lb(9)=0.95;
%  %lb(4)=T_eg;
%  ub=[E0,I0];
% % 
% % td=0:7:Ng*7;
% % l1=length(td)
% % NG-2
% % 
% %Q(tq)
% tl=0:42;
%  [paramfitL,resnormL] = lsqcurvefit(@Solve_sysL,paramguessL,tl,[E0,I0],lb,ub)
%  
%  E00=paramfitL(1)
%   I00=paramfitL(2)
%   
%   [tb,Solb] = ode45(@(t,X)EodeL(t,X,beta,phi,N),tl,[E00,I00]);
%   X0b=Solb(end,:)
%   
%  
%   Sb=Solb(:,1);
%  
%   figure(14)
%   plot(tb,Sb)
%   save('initval2b.mat','X0b');
 

% 
% 
% theta=thetafit;
% %theta=0.1;
% sigma=1/10;
% for i=1:NG-2
% R0G(i)=(p3*(CG(i+3)-(1-theta)*TG(i+3))+p2*(CG(i+2)-(1-theta)*TG(i+2))+p1*(CG(i+1)-(1-theta)*TG(i+1)))/(CG(i)-(1-theta)*TG(i));
% R0Gno(i)=(p3*(CG(i+3))+p2*(CG(i+2))+p1*(CG(i+1)))/(CG(i)-(1-theta)*TG(i));
% end

%for i=1:NS-2
% R0S(i)=(p3*(CS(i+3)-(1-theta)*TS(i+3))+p2*(CS(i+2)-(1-theta)*TS(i+2))+p1*(CS(i+1)-(1-theta)*TS(i+1)))/(CS(i)-(1-theta)*TS(i));
% R0Sno(i)=(p3*(CS(i+3))+p2*(CS(i+2))+p1*(CS(i+1)))/(CS(i)-(1-theta)*TS(i));
% end
% GavgR=mean(R0Gno);
% SavgR=mean(R0Sno)
% GavgRe=mean(R0G);

% t=1:length(TG);
% t2=1:length(TG)-3;
% ft = fittype( 'poly1' );
% fitresult = fit( t',qG', ft );
% fitresult1 = fit(t2',R0G,ft);
% fitresult2 = fit(t2',R0Gno,ft);
% figure (2)
% plot(1:Ng-3,R0G,'m-*',1:Ng-3,R0Gno,'r-*')
% hold on
% plot(fitresult1,'m--')
% hold on
% plot(fitresult2,'r--')
% 
% figure (3)
% plot(1:Ng,qG,'b-*',1:Ng,qGc,'g-*')
% hold on
% plot(fitresult,'b--')

% for i=1:length(ctran)
%     Td(i)=p*(fracCT)*((CTicdf(1:i)-CTicdf((1:i)-1))*Q(1:i)); 
%     while Td(i)>=C(i)
%       Td(i)=0.5*Td(i);
%      % jj=jj+1;
%     end  %traced infected from quarantined and traced to infected cdf%traced infected from quarantined and traced to infected cdf
% end
% 
% for i=length(ctran)+1:NQ
% %     CTicdf(ctran)-CTicdf(ctran-1)
%     Td(i)=p*(fracCT)*(CTicdf(ctran)-CTicdf(ctran-1))*Q(i-ctran); 
%    % jj=1;
% %    if Td(i)>=C(i)
% %      %uu=Td(i);
% %       Td(i)=C(i)*mean(Td(1:i-1)./C(1:i-1));
% %       
% %      % jj=jj+1;
% %     end
%     while Td(i)>=C(i)
%         uu=Td(i);
%       Td(i)=0.5*uu;
%      % jj=jj+1;
%     end%traced infected from quarantined and traced to infected cdf
% end

theta=0;
size(pmc)
size(C)
Td=pmc(1:length(C)).*C;
H=C-Td;

for i=1:Nr
    Re(i)=1/(H(i)+theta*Td(i))*((sicd(siran)-sicd(siran-1))*H(i+siran)+theta*(sicdCT(siran)-sicdCT(siran-1))*Td(i+siran)); 
    Reno(i)=1/(H(i)+theta*Td(i))*((sicd(siran)-sicd(siran-1))*H(i+siran)+(sicdCT(siran)-sicdCT(siran-1))*Td(i+siran));%traced infected from quarantined and traced to infected cdf
   % Renop(i)=1/(H(i))*((sicd(siran)-sicd(siran-1))*(H(i+siran).*(H(i+siran)+rho*Td(i+siran))./(H(i+siran)+(1-rho)*Td(i+siran))));
end

for i=1:Nc
    R0(i)=1/(C(i))*((sicd(siran)-sicd(siran-1))*C(i+siran)); 
end
% 
% figure(1)
% plot(1:Nr,Re,'b',1:Nr,Reno,'r',1:Nc,R0,'k')

 figure(12)
 plot(1:Nr,Re,'b-',1:Nr,Reno,'b--',1:Nr,Ret(1:Nr),'g-')


figure(2)
plot(1:Nd,C,'k-*',1:NQ,Td,'b-*',1:NQ,Qd,'c-*')


figure(6)
plot(1:Nd,C,'k-*',1:NQ,Td,'b-*',1:Nd,Idt,'k-',1:Nd,Idct,'b-')

 ItfigU(1:Nd,[C,Td,Idt,Iut+Idct,Idct])

end