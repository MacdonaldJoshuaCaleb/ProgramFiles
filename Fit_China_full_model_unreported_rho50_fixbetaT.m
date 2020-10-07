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
M(isnan(M))=0;
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
      u=0.75*conv(lognpdf(tv,1.57,0.65),logncdf(tv,0.74,0.64));
      u=u(1:length(tv))+0.25*logncdf(tv,1.57,0.65);
        
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

p=.01;  %prob. of infection/contact
theta=0.06;
% rho=0.5;
% for i=1:length(ctran)
%     Td(i)=p*(CTicdf(1:i)-CTicdf((1:i)-1))*Q(1:i); 
%     while Td(i)>=C(i)
%       Td(i)=0.5*Td(i);
%      % jj=jj+1;
%     end  %traced infected from quarantined and traced to infected cdf%traced infected from quarantined and traced to infected cdf
% end
% 
% for i=length(ctran)+1:NQ
% %     CTicdf(ctran)-CTicdf(ctran-1)
%     Td(i)=p*(CTicdf(ctran)-CTicdf(ctran-1))*Q(i-ctran); 
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
 N=1.6563e+08;
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
rho=0.75;
T=3;
T_e=3.71;
T_c=2.71;
Tg=T;
T_eg=T_e;
beta=4/N/Tg;
% Tg=4.64;
% T_eg=2.71;
% %phi1=mean(Td(5:14)./C(5:14));
% phi2=mean(Td(15:end)./C(15:end));
% phi1=0.5;
% sigma0=(1-p)/p
% c0=mean(Re(1:4))/(p*(1-phi1)*Tg+p*phi1*theta*Tg/T_eg*T_eg)/N;
% betag=c0*p;
% betaMg=c0*p*theta;
% Re1=(betag*(1-phi1)*Tg+be7aMg*phi1*T_eg)*N
alpha0=1/60;
nu0=0;
 Sc0=0;
 Ic0=0;
sqwt=0.05;
nu_c=0;
alpha_c=1/14;
fc=1;
function [dxdt]=Eode(t,X,beta,thetap,sigma,phi,nu,fc,thetac,p,alpha,rho)
        %beta=p(1);
        %T=p(2);
        
        betaM=thetap*beta;
        betac=thetac*beta;
        %T_c=0.5*T_e;
       % nu_c=0.5*nu;
        %
        %ps=1/sigma*(ddr-(1-p)/p*phi);
        %ps=0;
        %T_e=q(3);
      lambda=(beta*X(5)+betaM*X(6)+betac*X(11));
      
     dxdt(1)=-(1+(1-p)/p*phi*rho+sigma)*X(1)*(lambda)+alpha*X(2)+(1-fc)*alpha_c*X(9);
     dxdt(2)=(sigma)*X(1)*(lambda)-nu*X(2)*lambda-alpha*X(2)+fc*alpha_c*X(9);
     dxdt(3)= X(1)*(beta*(1-phi*rho)*X(5)+betaM*(1-phi*rho)*X(6)+betac*(1-phi)*X(11))- 1/tau*X(3);
     dxdt(4)= nu*(beta*(1-phi*rho)*X(5)+betaM*(1-phi*rho)*X(6)+betac*(1-phi)*X(11))*X(2)*lambda- 1/tau*X(4);
     dxdt(5)= 1/tau*X(3)-1/T*X(5);
     dxdt(6)=1/tau*X(4)-1/T_e*X(6);
     dxdt(7)=rho/T*X(5);
     dxdt(8)=rho/T_e*X(6);
     dxdt(9)=((1-p)/p*phi*rho)*X(1)*(lambda)-alpha_c*X(9)-nu_c*X(9)*X(5)*lambda;
     dxdt(10)= (nu_c*X(9)*lambda+(beta*(phi*rho)*X(5)+betaM*(phi*rho)*X(6)+betac*(phi)*X(11))*(X(1)+nu*X(2)))- 1/tau*X(10);
      dxdt(11)=1/tau*X(10)-1/T_c*X(11);
      dxdt(12)=1/T_c*X(11);
      dxdt(13)=1/T_c*X(11)-alpha_c*X(13);
     dxdt=dxdt';
end
% 
function[z] = Solve_sys(param,td)
        
      % pp=param(3);
       q=param(3:4);
       q2=param(5);
       nu=param(6);
       alpha=param(9);
       sigma=q2(1);
      % fc=q2(2);
       % beta=pp(1);
       % T=pp(2);
        thetac=param(7);
        p=param(8);
        phi=q(1);
        thetap=q(2);
        I0c=param(10);
       % T_e=q(3);
        I0cc=param(11);
        I0=param(1);
       % rho=param(13);
        Mc0=param(2);
%         
betaM=beta*thetap;    
betac=beta*thetac;
%         
% I0=C(1);
% Ic0=Td(1);
%Sc0=Q(1);
phim=1;
        E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1)+betaM*(1-phi*rho)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1)))*N;
        Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1)+betaM*(phi*rho)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1)))*N;
%         
%         
        S0=N;
%     
%      X0=[S0 Sc0 E0 Ec0 I0 0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0];
%      X0=[S0 Sc0 E0 Ec0 I0 0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0 0];
     X0=[S0 Sc0 E0 Ec0 I0 0 C(1) 0 (1-p)*Mc0 p*Mc0 I0c 0 I0cc];
%    
% td0=2:45;
     [t,Sol] = ode45(@(t,X)Eode(t,X,beta,thetap,sigma,phi,nu,fc,thetac,p,alpha,rho),td,X0);
     
%      figure(6)
%      plot(td0,Sol(:,5)+Sol(:,6),'r',td0,Sol(:,6),'b',1:Nd,C,'r*',1:NQ,Td,'b*')
%       
      CumC=(Sol(td,7)+Sol(td,8)+Sol(td,12));
      % CumSq=(Sol(tq,9)+Sol(tq,10)+Sol(tq,11))*sqwt;
       CumSq=(Sol(tq,9)+Sol(tq,10)+Sol(tq,11)+Sol(tq,13))*sqwt;
%       for j=2:length(td)
%           
%               Hh(j-1)=Sol(j,5)-Sol(j-1,5);
%               Tt(j-1)=Sol(j,7)-Sol(j-1,7);
%           
%       end
      z=[CumC',CumSq'];
end
% 
% T=4.64;
% %T=5
% T_e=3.71;
% T_c=2.71;
% Tg=T;
% T_eg=T_e;
T_cg=T_c;
%alpha=0;
%phi1=mean(Td(5:14)./C(5:14));
phi2=mean(Td(15:end)./C(15:end));
phi1=0.5;
 p=.05;
 p0=p;
sigma0=(1-p)/p+600;
c0=mean(Re(1:4))/(p*(1-phi1)*Tg+p*phi1*theta*Tg/T_eg*T_eg)/N;
betag=beta;
betaMg=beta*theta;
Re1=(betag*(1-phi1)*Tg+betaMg*phi1*T_eg)*N
rho0=.5;
 I0g=312/rho0;
 Mc0g=phi1/p*I0g;
 Mc0g=phi1*I0g;
 I0cg=p*Mc0g*.01;
 I0ccg=I0cg;

 thetac0=0.05;
ps0=1/sigma0*(ddr-(1-p0)/p0*phi1)
 paramguess=[I0g,Mc0g,phi1,theta,sigma0,nu0,thetac0,p0,alpha0,I0cg,I0ccg];
 lb=zeros(1,length(paramguess));
 %lb(3)=betag*.25;
 %lb(end-4)=0.95;
 %lb(4)=.2;
 ub=[10*I0g,10*Mc0g,1,.15,4000,1,.1,.1,.1,10*I0cg,10*I0ccg];
% 
% td=0:7:Ng*7;
% l1=length(td)
% NG-2
% 
%Q(tq)
 [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[CC,M(tq)*sqwt],lb,ub)
% [I0g,Mc0g,betag,phi1,theta,sigma0,1/2,nu0,thetac0,p0];
save('parametersNp2FLwinc2ccR2e.mat','paramfit');
% betaf=paramfit(3)
betaf=beta
 phif=paramfit(3)
  %fc=paramfit(7)
  nuf=paramfit(6)
  alpha=paramfit(9);
  alphaf=alpha
 Tfit=T
  Tefit=T_e
  Rfit=betaf*Tfit*N
  sigmaf=paramfit(5)
 % ps=1/sigmaf*(ddr-(1-p)/p*phif)
 betacf=paramfit(7)*betaf
  thetac=paramfit(7)
 p=paramfit(8)
% Rmfit=paramfit(6)*paramfit(7)
 thetafit=paramfit(4)
 fracCT=1
 I0=paramfit(1)
Mc0=paramfit(2)
I0c=paramfit(10)
       % T_e=q(3);
        I0cc=paramfit(11)
 beta=betaf;
 phi=phif;
psi=sigmaf+(1-p)/p*phi
 %T_e=Tefit;
% rho=paramfit(13)
 T=Tfit;
betaM=beta*thetafit; 
nu=nuf;
S0=N;
phim=1;
 E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1)+betaM*(1-phi*rho)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(1-phim)*T*I0/(1/tau*T+1)+betaM*(1-phim)*T_e*Ic0/(1/tau*T_e+1)))*N;
        Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1)+betaM*(phi*rho)*T_e*Ic0/(1/tau*T_e+1)+nu*(beta*(phim)*T*I0/(1/tau*T+1)+betaM*(phim)*T_e*Ic0/(1/tau*T_e+1)))*N;
%                 
%     X0=[S0 Sc0 E0 Ec0 I0 Ic0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0]
%      X0=[S0 Sc0 E0 Ec0 I0 0 C(1) Td(1) (1-p)*Mc0 p*Mc0 0 0 0];
 X0=[S0 Sc0 E0 Ec0 I0 0 C(1) 0 (1-p)*Mc0 p*Mc0 I0c 0 I0cc];
        
 Rmfit=betacf*T_c*N
 Refit=(1-phif)*Rfit+phif*Rmfit
% betaf-betag
% phif-phig
% 
% betaf-betag
% phif-phig
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
 
 [tn,Soln] = ode45(@(t,X)Eode(t,X,beta,thetafit,sigmaf,phi,nu,fc,thetac,p,alpha,rho),td,X0);
 
 figure(10)
 plot(tn,Soln(:,5)+Soln(:,6)+Soln(:,11))
  figure(11)
 %plot(tn,Soln(:,9)+Soln(:,10)+Soln(:,11))
plot(tn,Soln(:,9)+Soln(:,10)+Soln(:,11)+Soln(:,13))
 
  St=Soln(:,1);
  Sct=Soln(:,9);
  Sqt=Soln(:,2);
 It=Soln(:,5);
 Ict=Soln(:,11);
  Idt=(diff(Soln(:,7))+diff(Soln(:,8))+diff(Soln(:,12)));
 Idct=diff(Soln(:,12));
 Itt=It+Ict;
 
 pm=Ict./ (It+Ict);
  pmc=diff(Soln(:,12))./(diff(Soln(:,7))+diff(Soln(:,8))+diff(Soln(:,12)));
  
  Ret=zeros(length(tn),1);
  
  for i=1:length(tn)
 F=[0, 0,0, (1-phi)*beta*St(i),(1-phi)*betac*St(i),(1-phi)*betaM*St(i); 0,0,0,phi*beta*(St(i)+nu*Sqt(i)),phi*betac*(St(i)+nu*Sqt(i)),phi*betaM*(St(i)+nu*Sqt(i)); 0,0,0, (1-phi)*beta*nu*Sqt(i),(1-phi)*betac*nu*Sqt(i),(1-phi)*nu*betaM*Sqt(i); 0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0];
V=[1/tau,0,0,0,0,0;0,1/tau,0,0,0,0;0,0,1/tau,0,0,0;-1/tau,0,0,1/T,0,0;0,-1/tau,0,0,1/T_c,0;0,0,-1/tau,0,0,1/T_e];

         M=F*inv(V);
       Ret(i)=max(eig(M));
  end
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

theta=Rmfit/Rfit
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

 figure(1)
 plot(1:Nr,Re,'b-',1:Nr,Reno,'b--',1:Nr,Ret(1:Nr),'g-')


figure(2)
plot(1:Nd,C,'k-*',1:NQ,Td,'b-*',1:NQ,Qd,'c-*')


figure(6)
plot(1:Nd,C,'k-*',1:NQ,Td,'b-*',1:Nd,Idt,'k-',1:Nd,Idct,'b-')


end