ags = linspace(.03,4,400);
bin = zeros(length(ags),5);
for w = 1:length(ags)
   bin(w,:) = FitSimpleModelChinaDailyGammaQLoopfunc(ags(w));
   fprintf('%i of %i\n',w,length(ags))
end

I0s = bin(:,1);
R0bs = bin(:,2);
sigmas = bin(:,3);
phis = bin(:,4);

Resids = bin(:,5);

figure
plot(ags,I0s,'linewidth',2)
title('I_0 as fxn of ag, China')



figure
plot(ags,R0bs,'linewidth',2)
title('R0_b as fxn of ag, China')


figure
plot(ags,sigmas,'linewidth',2)
title('\sigma as fxn of ag, China')

figure
plot(ags,phis,'linewidth',2)
title('\phi as fxn of ag, China')


figure
plot(ags,Resids,'linewidth',2)
title('Residual as fxn of ag')

%%
close all
load('tq.mat')
load('paramfitChinaGammaQ.mat')
load('paramfitChinaWeibullQ.mat')
%tq = linspace(.1,48);
ag = paramfitChinaGammaQ(5)*14
bg = 14/ag
bw = paramfitChinaWeibullQ(5)
aw = 14/gamma(1+(1/bw))

figure
plot(tq,gamcdf(tq,ag,bg),tq,wblcdf(tq,aw,bw),tq,expcdf(tq,(14)),'linewidth',2)
legend('gamma','weibull','exponential')
title('CDF')


figure
plot(tq,gampdf(tq,ag,bg),tq,wblpdf(tq,aw,bw),tq,exppdf(tq,(14)),'linewidth',2)
legend('gamma','weibull','exponential')
%ylim([0, .5])
title('PDF')