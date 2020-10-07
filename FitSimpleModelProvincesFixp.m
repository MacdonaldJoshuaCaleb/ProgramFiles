function []= FitSimpleModelProvincesFixp
close all
warning off;
% stop/start figures from displaying in matlab, switch 'off' to 'on to display
% import data 
ChinaData = readtable('ChinaDataReshaped.csv');
ChinaData = table2array(ChinaData(:,2:end));
ChinaData(:,30) = [];
Province = readtable('ProvinceNames.csv');
Province = table2array(Province);
Province = string(Province);
Province(30) = [];
ProvincePops = readtable('ChinaProvincePops.csv');
ProvincePops = table2array(ProvincePops);
ProvincePops(30) = [];
% national total cases
CCt = ChinaData(:,end-2);
% national total quarintined 
QCt = ChinaData(:,end);
% national current quarintined
%Mt = table2array(ChinaData(:,end));
%Qtd = Mt(2:end);
% national daily cases
Ctd = diff(CCt);
%tdt=1:length(CCt);
% storage bin for fit parameteres 
FitparamsProvFixp = zeros(length(ProvincePops),4);

T=4.64;
T_e=2.71;
tau=3;
rho=1;
sqwt=.01;
alphac=1/14;
p = .06;
Mc0 = 0.00172837295386061;
%beta = 6/T;
% ODE to be fit
%function [dxdt]=Eode(t,X,beta,sigma,phi,p,rho,T,N)
%    psi = sigma + ((1-p)/p)*phi*rho;
%    dxdt(1)=-(1+psi)*X(1)*((beta/N)*X(3));
%    dxdt(2)= (beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2);
%    dxdt(4)= (beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4);
%    dxdt(3)= 1/tau*X(2)-1/T*X(3);
%    dxdt(5)=1/tau*X(4)-1/T_e*X(5);
%    dxdt(6)=(beta/N)*(phi*rho)*X(1)*X(3);
%    dxdt(7)= rho/T*X(3);
%    dxdt(8)=1/T_e*X(5);
%    dxdt(9)=(beta/N)*(phi*rho)/p*(1-p)*X(1)*X(3)-alphac*X(9)+1/T_e*X(5);
%    dxdt=dxdt';
%end

function X0 = IC(param,QC,C,M,tq,I0)
        I0 = param(1);
        beta = param(2);
        sigma = param(3);
        phi = param(4);
        %p=param(4);
        N = param(5);
        S0=N;
        %E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1));
        E0=beta*(1-phi*rho)*tau*I0;
        %Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1));
        Ec0=(beta*(phi*rho)*tau*I0);
        Ic0=phi*rho*I0;
        M0 = M(tq(1));
        X0=[S0 E0 I0 p*tau*M0/(tau+1) p*M0/(tau+1) p*M0 (1-phi)*C(1) phi*rho*C(1) M0];

        

        %X0=[S0 E0 I0 Ec0 Ic0 QC(1) C(1) Td(1) Mc0];        
end

% Objective Function 
function[z] = Solve_sys(param,td)
        if its == 10
            sqwt = .2;
        end
        if its == 12
            sqwt = .2;
        end
        if its == 14
            sqwt = .05;
        end
        if its == 20
            sqwt = .1;
        end
        if its == 23
            sqwt = .05;
        end
        I0 = param(1);
        beta = param(2);
        sigma = param(3);
        phi = param(4);
        %p=param(4);
        S0 = param(5);
        N = S0;
        

      %  E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1));
       %  E0=beta*(1-phi*rho)*tau*I0;
       % Ec0=(beta*(phi*rho)*tau*I0);
       % Ic0=phi*rho*I0;
      
        %S0=N;
        X0 = IC(param,QC,C,M,tq,I0g);
        S0 = X0(1);
        %N = S0;
        psi = sigma + ((1-p)/p)*phi*rho;
        Eode = @(t,X) [-(1+psi)*X(1)*((beta/N)*X(3)),(beta/N)*(1-phi*rho)*X(1)*X(3)- 1/tau*X(2),...
                       1/tau*X(2)-1/T*X(3),(beta/N)*(phi*rho)*X(1)*X(3)- 1/tau*X(4),...
                       1/tau*X(4)-1/T_e*X(5),(beta/N)*(phi*rho)*X(1)*X(3),...
                       rho/T*X(3),1/T_e*X(5),(beta/N)*(phi*rho)/p*X(1)*X(3)-alphac*X(9)]';

       [t,Sol] = ode45(@(t,X)Eode(t,X),td,X0);
       
   

      Cdd(2:end)=Sol(2:end,7)+Sol(2:end,8)-(Sol(1:end-1,7)+Sol(1:end-1,8));
      
        CumSq= Sol(tq-(tq(1)-1),9)*sqwt;
    

          

        CumSq(1) = M(tq(1))*sqwt;


      Cds=Cdd(tds);
      z=[Cds',CumSq'];
      
      
end



% [I0, Mc0,psi,phi,p] from previous fitting
% we will fix all parameters except Beta, phi, sigma, I0
% psi = sigma + ((1-p)/p)*phi
%paramfit1 = [778.050470654426,0.00172837295386061,1244.68246300590,0.376916586571512,0.0816546156482878];
%p = paramfit1(5);

for its = 1:length(Province)
    set(0,'DefaultFigureVisible','on')
        if its == 10
            sqwt = .2;
        end
        if its == 12
            sqwt = .2;
        end
        if its == 14
            sqwt = .05;
        end
        
        if its == 20
            sqwt = .1;
        end
        if its == 23
            sqwt = .05;
        end

    % retrieve province population
    N = ProvincePops(its);
    % retrieve province name
    Name = Province(its);
    % provincial total cases
    CC = ChinaData(:,its);
    CC = CC(CC > 0);
    td=1:length(CC);
    zers = length(CCt)-length(CC)+1;
    % provinicial total quarintined estimate 
    %QC = ceil((CC/CCt(zers:end))*QCt(zers:end));
    f = @(t) exp(-(1/14).*t);
    
    % provincial current quarintined estimate
    %Qd = diff(QC);
    %M = [QC(1);Qd'];
    % provincial daily cases
    Cd = diff(CC);
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
if its == 14 
ind0=find(Cd==0);
if ind0(end) == length(Cd)
    ind0 = ind0(1:end-1);
end

   if isempty(ind0)==0 
    Cd(ind0)=1/3*Cd(ind0-1)+1/3*Cd(ind0+1);
    Cd(ind0-1)=2/3*Cd(ind0-1);
    Cd(ind0+1)=2/3*Cd(ind0+1);
   end
end

if its == 23 
    
ind0=find(Cd==0);
if ind0(end) == length(Cd)
    ind0 = ind0(1:end-1);
end

   if isempty(ind0)==0 
    Cd(ind0)=1/3*Cd(ind0-1)+1/3*Cd(ind0+1);
    Cd(ind0-1)=2/3*Cd(ind0-1);
    Cd(ind0+1)=2/3*Cd(ind0+1);
   end
end

while Cd(1) >= Cd(2) || Cd(1) >= Cd(3)
    Cd(2) = Cd(2) + Cd(1)/3;
    Cd(3) = Cd(2) + Cd(1)/3;
    Cd(1) = Cd(1)/3;
end
if its == 12
    temp = Cd(49);
    Cd(49) = 0;
    Cd(37) = temp;
end
    

%if ind0(end) == length(Cd)
%    if isempty(ind0(1:end-1))==0
%        Cd(ind0(1:end-1))=1/3*Cd(ind0(1:end-1)-1)+1/3*Cd(ind0(1:end-1)+1);
%        Cd(ind0(1:end-1)-1)=2/3*Cd(ind0(1:end-1)-1);
%        Cd(ind0(1:end-1)+1)=2/3*Cd(ind0(1:end-1)+1);
%    end
%    if isempty(ind0(end)) == 0
%        Cd(ind0(end))= 1/3*Cd(ind0(end)-1)+1/3*Cd(ind0(end)-2);
%    end
%end

    
  
    
    v1 = zeros(1,length(td));
    v2 = zeros(1,length(td));
    v1(1) = Cd(1);
    v2(1) = Ctd(1);
    Mt=ChinaData(:,end-1);
    QC(1) = (v1(1)/v2(1))*Mt(zers);
    %length(f(td(2:end).*Ctd(zers+1:end)')
    
    v1(2:end) = (f(td(2:end)).*Cd');
    v2(2:end) = (f(td(2:end)).*Ctd(zers:end)');



    for w = 2:length(v1)
        QC(w) = trapz(v1(1:w))/trapz(v2(1:w))*Mt(zers+w-1);
    end
    %if its == 23
    %    QC(1) = QC(4)-30;
    %    QC(2) = QC(4) - 20;
    %    QC(3) = QC(4) - 10;
    %end
    
    M = QC;


    %M = QC';
    tq = find(QC > 0);
    if its == 23
        M((tq(4))) = 2*M(tq(5))-M(tq(6));
        M((tq(3))) = 2*M(tq(4))-M(tq(5));
        M((tq(2))) = 2*M(tq(3))-M(tq(4));
        M((tq(1))) = 2*M(tq(2))-M(tq(3));
    end
    

    M(tq(1:4))
        
    C = Cd(1:end);
    [mx,im]=max(Cd); 
    % account for discrepancy in data for Hubei Province
    if its ~= 14
        tds = td(1:end-1);
    end
    if its == 14
        tds=td([1:im-2,im+2:length(Cd)]);
        C = C(tds);
        
    end
    if its == 32
        tds = td([1:29,31:length(Cd)]);
        C = C(tds);
    end
    Cdd = CC;
    Cdd(1)=Cd(1);
    Td=.1609*C;
    I0g = (C(1)/rho)*T;
    % Hong Kong data very noisy, come back to this
    %if its == 13
    %    I0g = 0;
    %end
    % fit values from entire country fit with exception of I0
    %paramguess = [I0g,	(1.66014398889413e-08)*(1e8),	1407.12998368299,	0.391575667314593,	0.0763792803566279,N];
    paramguess = [I0g,(1.66014398889413e-08)*(1e8),	1407.12998368299,	0.391575667314593,N];
    lb = [I0g*.5,4/T,1.1500e+03,.01,N-eps];
    %if its == 13
    %    lb(1) = 0;
    %end
    ub = [2*I0g,6/T,inf,1,N+eps];
    %paramguess = (ub-lb)*.5;
    % reshape data to get into rows for fitting 
    %M = M;
    C = C';
    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramguess,td,[C,M(tq)*sqwt],lb,ub);
    lb(1) = paramfit(1)*.5;
    ub(1) = paramfit(1)*2;
    [paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C,M(tq)*sqwt],lb,ub);
    %[paramfit,resnorm] = lsqcurvefit(@Solve_sys,paramfit,td,[C,M(tq)*sqwt],lb,ub);
    % store fit parameters 
    FitparamsProvFixp(its,1:4) = paramfit(1:4);
    
    % post process plot generation 
    I0 = paramfit(1)
    beta = paramfit(2)
    sigma = paramfit(3)
    phi = paramfit(4)
    p
    %psi = sigma + ((1-p)/p)*phi*rho
    %FitparamsProv(its,6) = psi;
    S0=N;
    %Mc0=paramfit1(2);

    %E0=(beta*(1-phi*rho)*T*I0/(1/tau*T+1));
    %Ec0=(beta*(phi*rho)*T*I0/(1/tau*T+1));
    %Ic0=phi*rho*I0;
    %X0=[S0 E0 I0 Ec0 Ic0 p*QC(1) C(1) phi*rho*C(1) (1-p)*M(tq(1))];
    zFit=Solve_sys(paramfit,td);
    CFit=zFit(1:length(tds));
    QFit=zFit(length(tds)+1:end);
    %QFits(its,:) = QFit;
    set(0,'DefaultFigureVisible','off')
    figure 
    plot(tds,CFit,'b-',td(1:end-1),Cd,'r*','linewidth',2)
    str = strcat('New Reported Cases,',{' '},Province(its),{' '},'Province');
    title(str)
    xlabel('Days')
    ylabel('Count')
    ylim([0 inf])
    baseFileName = sprintf('Cases%d.png', its);
    fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\ProvinceFittingPlots';
    saveas(gca, fullfile(fname, baseFileName), 'png');
    figure 
    plot(tq,QFit/(sqwt),'b-',tq,M(tq),'r*','linewidth',2)
    str = strcat('Estimated Quarintined Individuals,',{' '},Province(its),{' '},'Province');
    title(str)
    xlabel('Days')
    ylabel('Count')
    ylim([0 inf])
    baseFileName = sprintf('Quars%d.png', its);
    fname = 'C:\Users\macdo\Dropbox\COVID-19_modelfitting_data\ProvinceFittingPlots';
    saveas(gca, fullfile(fname, baseFileName), 'png');
    
    
    
    % clear reused arrays
    clear CC;
    clear QC;
    clear Qd;
    clear M;
    clear Cd;
    clear td;
    clear tq;
    clear tds;
    clear C;
    clear Cdd;
    clear Td;
    clear N;
    clear t;
    clear Sol;
    clear v1;
    clear v2;
    clear tc;
    clear ind0;
    % iteration tracker
    fprintf('fit % i of % i\n',its,length(Province))
end
% save fit parameters 
FitparamsProvFixp(:,2) = FitparamsProvFixp(:,2)*T;
FitparamsProvFixp(13,:) = [];
save('FitparamsProvFixp','FitparamsProvFixp')
%save('QFits.mat','QFits')
end