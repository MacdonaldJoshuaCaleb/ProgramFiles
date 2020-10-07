rho = 1;
T=4.64;
load('paramfitChina.mat')
paramfit = paramfitChina;

load('FitparamsSynthChina.mat')
temp = FitparamsSynth;
temp(:,2)= temp(:,2)*T;
tempS = sort(temp);
tempS = tempS(251:10000-251,:);
CILB = zeros(1,4);
CIUB = zeros(1,4);
Xbar = zeros(1,4);
S = zeros(1,4);
ARE = zeros(1,4);
MLE = paramfit;


for j = 1:4
    CILB(j) = tempS(1,j);
    CIUB(j) = tempS(end,j);
    Xbar(j) = mean(temp(:,j));
    S(j) = std(temp(:,j));
    ARE(j) = sum(abs(temp(:,j)-MLE(j))/MLE(j))*(1/10000)*(100);
end
CILB = CILB';
CIUB = CIUB';
Xbar = Xbar';
S = S';
ARE = ARE';
MLE = MLE';
T2 = [CILB,CIUB,MLE,Xbar,S,ARE];
latex_table = latex(vpa(sym(T2),5))