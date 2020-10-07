function []= FitSimpleModel
warning off;
filename = 'data-preview_cn_apr.xlsx';
xlRange = 'E34:BJ37';
dataArray = xlsread(filename,xlRange);
CC=dataArray(1,:);
QC=dataArray(2,:);
M=dataArray(3,:);
Qd=M(2:end);
Cd=diff(CC);
td=1:length(CC);
Qd(isnan(Qd))=0;
QC(isnan(QC))=0;
tq=find(QC>0);
C=Cd(1:end);
C=C';
Td=.1609*C;
T=4.64;
T_e=2.71;
N=1e+08;
tau=3;
rho=1;
sqwt=0.05;
alphac=1/14;

% [I0, Mc0,psi,phi,p] from previous fitting
% we will fix all parameters except Beta, phi, sigma, I0
% psi = sigma + ((1-p)/p)*phi
paramfit1 = [778.050470654426,0.00172837295386061,1244.68246300590,0.376916586571512,0.0816546156482878];
p0 = paramfit1(5);
phi0 = paramfit1(4);


 function [dxdt]=Eode(t,X,beta,sigma,phi,p,rho,T)
    psi = sigma + ((1-p)/p)*phi*rho;
    dxdt(1)=-(1+psi)*X(1)*(beta*X(3));
    dxdt(2)= beta*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2);
    dxdt(4)= beta*(phi*rho)*X(1)*X(3)- 1/tau*X(4);
    dxdt(3)= 1/tau*X(2)-1/T*X(3);
    dxdt(5)=1/tau*X(4)-1/T_e*X(5);
    dxdt(6)=beta*(phi*rho)*X(1)*X(3);
    dxdt(7)= rho/T*X(3);
    dxdt(8)=1/T_e*X(5);
    dxdt(9)=beta*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5);
    dxdt=dxdt';
 end

 
        
     
     


% 
function[z] = Solve_sys(param,td,beta)
        I0 = param(1);
        %beta = param(2);
        sigma = param(2);
        phi = param(3);
        p=param(4);
        Mc0=paramfit1(2);
        S0=N;

        E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N;
        Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
        Ic0=phi*rho*I0;
        S0=N;
        X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
 
       [t,Sol] = ode45(@(t,X)Eode(t,X,beta,sigma,phi,p,rho,T),td,X0);

      CumC=Sol(td,7)+Sol(td,8);
      CumSq=(Sol(tq,4)+Sol(tq,5)+Sol(tq,9))*sqwt;

      z=[CumC',CumSq'];
end

R0bvec=linspace(4,12,41);
R0vec=zeros(size(R0bvec));
xvec=zeros(size(R0bvec));
phivec=zeros(size(R0bvec));
Residvec=zeros(size(R0bvec));

I0g=paramfit1(1)/rho;
sigmag = paramfit1(3) - ((1-p0)/p0)*phi0*rho;
nu=0;
for i=1:length(R0bvec)
beta=R0bvec(i)/N/T;
f=@(param,td)Solve_sys(param,td,beta);
pg=R0bvec(i)/6*p0;
 phig = R0bvec(i)/6*phi0;
%  pg=p0;
%  phig = phi0;
%phig=0.05;
%p=.01;


paramguess=[I0g,sigmag,phig,pg];
 lb=[0,0,.05,0];
 ub = [10*I0g,2000-((1-p0)/p0)*phi0,1,1];

 [paramfit,resnorm] = lsqcurvefit(f,paramguess,td,[CC,M(tq)*sqwt],lb,ub)
% 
%save('fit_paramsp.mat','paramfit');
 
phivec(i)=  paramfit(3);
Residvec(i)=resnorm;
R0vec(i)=R0bvec(i)*(1-phivec(i));
Re=R0vec(i);
I0=paramfit(1);

sigma = paramfit(2)
phi = paramfit(3)
p=paramfit(4)
S0=N;
psi=sigma+(1-p)/p*phi;


E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1))*N;
  fun=@(x)  log(x)-(1+psi)*(Re*(x-1+(psi/(-nu+1+psi))*((x)^(nu/(1+psi))-x))-(beta*T*(E0+I0)));
         %fun=@(x)  log(x)-(1+psi)*Re*(x-1);
        % Xs=fzero(fun,x0);
         Xs=bisection(fun,0,1);
          Cs=(1-psi/(-nu+1+psi)*((-nu+1)/psi*Xs+Xs^(nu/(1+psi))))*S0;  %Xs=1-Sinf/St
         Cs=Cs+C(1)*((1-phi)/rho+phi);
         
         xvec(i)=Cs;

end


% Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1))*N;
% Ic0=phi*rho*I0;
% X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];
% 
% 
%  zFit=Solve_sys(paramfit,td);
% % 
%  CFit=zFit(1:length(td));
%  DCFit = zeros(length(CFit));
%  DCFit(1) = 0;
%  for j = 2:length(CFit)
%      DCFit(j) = CFit(j) - CFit(j-1);
%  end
%  DCFit = DCFit(2:length(td));
%  QFit=zFit(length(td)+1:end);
% % 
 figure(1)
 plot(R0bvec,phivec,'b-')
 
 
 figure(2)
 plot(R0bvec,Residvec,'r*')
 
 figure(3)
 plot(R0bvec,R0vec,'k-')
 
  figure(4)
 plot(R0bvec,xvec,'k-')

%  
%  figure 
%  plot(tq,QFit/sqwt,'b-',td,M,'r*')
%  title('quarantined')
%  
% figure
% plot( td(2:end),DCFit,'b',td(2:end),C,'r*')
% 
% R0b = beta*T*N
% R0 = R0b*(1 - phi*rho)
% repos = [R0,R0b];
% save('repo_nums_p.mat','repos');
%  

 end

