Ts = linspace(4.13,5.1,40);
Resid = zeros(length(Ts),1);
for w = 1:length(Ts)
    Resid(w) = FitSimpleModelSerialIntVar(Ts(w));
end
set(0,'DefaultFigureVisible','on')
figure
plot(Ts,Resid,'linewidth',2)
ylabel('Residual')
xlabel('T')