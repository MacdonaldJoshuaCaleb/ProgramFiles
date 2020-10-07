tq = linspace(.1,14);
close all
figure
plot(tq,gamcdf(tq,2,4.64/2),tq,expcdf(tq,4.64),tq,logncdf(tq,1.23,.79),'linewidth',2)
title('Considered CDFs for T')
legend('Erlang','exponential','lognormal')
figure
plot(tq,gampdf(tq,2,4.64/2),tq,exppdf(tq,4.64),tq,lognpdf(tq,1.23,.79),'linewidth',2)
title('Considered PDFs for T')
legend('Erlang','exponential','lognormal')

figure
plot(tq,gamcdf(tq,2,2.71/2),tq,expcdf(tq,2.71),tq,logncdf(tq,.77,.67),'linewidth',2)
title('Considered CDFs for T_e')
legend('Erlang','exponential','lognormal')
figure
plot(tq,gampdf(tq,2,2.71/2),tq,exppdf(tq,2.71),tq,lognpdf(tq,.77,.67),'linewidth',2)
title('Considered PDFs for T_e')
legend('Erlang','exponential','lognormal')

figure
plot(tq,gamcdf(tq,2,3/2),tq,expcdf(tq,3),'linewidth',2)
title('Considered CDFs for \tau')
legend('Erlang','exponential')
figure
plot(tq,gampdf(tq,2,3/2),tq,exppdf(tq,3),'linewidth',2)
title('Considered PDFs for \tau')
legend('Erlang','exponential')

%%