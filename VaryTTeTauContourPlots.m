Ts = linspace(4.13,5.1,20);
T_es = linspace(2.08,3.31,20);
taus = [2,3,4];
[X,Y,Z] = meshgrid(Ts,T_es,taus);
Fits = zeros(20,20,3,13);
for w = 1:20
    for s = 1:20
        for t = 1:3
            Fits(w,s,t,:) = FitSimpleModelVar(X(w,s,t),Y(w,s,t),Z(w,s,t));
            if Fits(w,s,t,2) > 6
                Fits(w,s,t,2) = 6;
            end
            if Fits(w,s,t,6) > 6
                Fits(w,s,t,6) = 6;
            end
        end
    end
    fprintf('FITTING %i OF %i\n',(w*20*3),20*20*3)
end



%%
close all

figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,1),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,1),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,1),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit I_0 for China Less Hubei')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,2));
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,2));
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,2),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit R_{0,b} for China Less Hubei')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,3),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,3),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,3),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit \sigma for China Less Hubei')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,4),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,4),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,4),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit \phi for China Less Hubei')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,7),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,7),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,7),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit \sigma for Hubei Province')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,6));
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,6));
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,6),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit R_{0,b} for Hubei Province')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,5),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,5),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,5),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit I_0 for Hubei Province')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,8),5);
title('\tau=4')
%colormap(ax1,winter)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,8),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,winter)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,8),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,winter)
colorbar
sgtitle('Fit \phi for Hubei Province')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,9),5);
title('\tau=4')
%colormap(ax1,sum)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,9),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,summer)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,9),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,summer)
colorbar
sgtitle('Residual')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,10),5);
title('\tau=4')
%colormap(ax1,sum)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,10),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,summer)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,10),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,summer)
colorbar
sgtitle('Final Outbreak Size for China Less Hubei')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,11),5);
title('\tau=4')
%colormap(ax1,sum)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,11),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,summer)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,11),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,summer)
colorbar
sgtitle('Final Outbreak Size Reduction for China Less Hubei (%)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,13),5);
title('\tau=4')
%colormap(ax1,sum)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,13),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,summer)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,13),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,summer)
colorbar
sgtitle('Final Outbreak Size Reduction for Hubei Province (%)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
tiledlayout(3,1)
ax1 = nexttile;
contourf(X(:,:,3),Y(:,:,3),Fits(:,:,3,12),5);
title('\tau=4')
%colormap(ax1,sum)
colorbar
ax2 = nexttile; 
hold on
contourf(X(:,:,2),Y(:,:,2),Fits(:,:,2,12),5);
plot(4.64,2.71,'k*','linewidth',2)
hold off
title('\tau=3')
ylabel('T_e')
%colormap(ax2,summer)
colorbar
ax3 = nexttile;
contourf(X(:,:,1),Y(:,:,1),Fits(:,:,1,12),5);
title('\tau=2')
xlabel('T')
%colormap(ax3,summer)
colorbar
sgtitle('Final Outbreak Size for Hubei Province')


