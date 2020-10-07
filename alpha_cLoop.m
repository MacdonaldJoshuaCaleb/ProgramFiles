close all
alphacs = linspace(1/16,1/4,50);
bin = zeros(length(alphacs),13);
for w = 1:length(alphacs)
   bin(w,:) = FitSimpleModelVarQuar(alphacs(w));
   fprintf('ITERATION %i OF %i\n',w,length(alphacs))
end

CLHI0s = bin(:,1);
CLHR0bs = bin(:,2);
CLHsigmas = bin(:,3);
CLHphis = bin(:,4);

HI0s = bin(:,5);
HR0bs = bin(:,6);
Hsigmas = bin(:,7);
Hphis = bin(:,8);

Resids = bin(:,9);
CLHFS=bin(:,10);
CLHFSred=bin(:,11);
HFS=bin(:,12);
HFSred=bin(:,13);

% I0
figure
plot(alphacs,CLHI0s,'b','linewidth',2)
title('China Less Hubei')
ylabel('fit values for I_0')
xlabel('\alpha_c')
xlim([1/16,1/4])
figure
plot(alphacs,HI0s,'b','linewidth',2)
title('Hubei')
ylabel('fit values for I_0')
xlabel('\alpha_c')
xlim([1/16,1/4])
% R0
figure
plot(alphacs,CLHR0bs,'b','linewidth',2)
title('China Less Hubei')
ylabel('fit values for R_{0,b}')
xlabel('\alpha_c')
xlim([1/16,1/4])
ylim([.99999*min(CLHR0bs),1.00001*max(CLHR0bs)])
figure
plot(alphacs,HR0bs,'b','linewidth',2)
title('Hubei')
ylabel('fit values for R_{0,b}')
xlabel('\alpha_c')
xlim([1/16,1/4])
ylim([.99999*min(HR0bs),1.00001*max(HR0bs)])
% sigma
figure
plot(alphacs,CLHsigmas,'b','linewidth',2)
title('China Less Hubei')
ylabel('fit values for \sigma')
xlabel('\alpha_c')
xlim([1/16,1/4])

figure
plot(alphacs,Hsigmas,'b','linewidth',2)
title('Hubei')
ylabel('fit values for \sigma')
xlabel('\alpha_c')
xlim([1/16,1/4])
% phi

figure
plot(alphacs,CLHphis,'b','linewidth',2)
title('China Less Hubei')
ylabel('fit values for \phi')
xlabel('\alpha_c')
xlim([1/16,1/4])

figure
plot(alphacs,Hphis,'b','linewidth',2)
title('Hubei')
ylabel('fit values for \phi')
xlabel('\alpha_c')
xlim([1/16,1/4])


% residual 

figure
plot(alphacs,Resids,'b','linewidth',2)
ylabel('Residual')
xlabel('\alpha_c')
xlim([1/16,1/4])

% Outbreak Size
figure
plot(alphacs,CLHFS,'b','linewidth',2)
title('China Less Hubei')
ylabel('Final Outbreak Size')
xlabel('\alpha_c')
xlim([1/16,1/4])

figure
plot(alphacs,HFS,'b','linewidth',2)
title('Final Size as fxn of R_{0,b}, Hubei')
title('Hubei')
ylabel('Final Outbreak Size')
xlabel('\alpha_c')
xlim([1/16,1/4])

figure
plot(alphacs,CLHFSred,'b','linewidth',2)
title('China Less Hubei')
ylabel('Final Outbreak Size Reduction (%)')
xlabel('\alpha_c')
xlim([1/16,1/4])

figure
plot(alphacs,HFSred,'b','linewidth',2)
title('Hubei')
ylabel('Final Outbreak Size Reduction (%)')
xlabel('\alpha_c')
xlim([1/16,1/4])
