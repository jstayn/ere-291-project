data = csvread('sensitivity-output.csv');

x = data(:,1);
y = data(:,2);
c = data(:,3);
r = data(:,4);

greenX = x(r==1);
greenY = y(r==1);
greenC = c(r==1);

redX = x(r==3);
redY = y(r==3);
redC = c(r==3);

blackX = x(r==4);
blackY = y(r==4);
blackC = c(r==4);

logicBC = (0 < blackC) & (blackC < 50);
blackX1 = blackX(logicBC);
blackY1 = blackY(logicBC);
blackC1 = blackC(logicBC);

figure(1)
scatter3(greenX',greenY',greenC',10,'g','filled')
hold all
scatter3(redX',redY',redC',5,'r','filled')
scatter3(blackX1',blackY1',blackC1',5,'black','filled')
xlabel('CMH [m^3 / hr]')
ylabel('CO_2 [ppm]')
zlabel('OpEx Cost [$]')
title('Starting Point Sensitivity Analysis, Day One OpEx')
legend('Optimal','Infeasible','Time-out')

%{
    figure(2)
    scatter3(greenX',greenY',greenR',3,'g','.')
    hold all
    scatter3(redX',redY',redR',3,'r','.')
    scatter3(blackX',blackY',blackR',3,'black','.')
    xlabel('Tank Radius [m]')
    ylabel('Number of Tanks')
    zlabel('Optimal Radius [m]')
    title('Starting Point Sensitivity Analysis, Tank Radius')
    legend('Optimal','Infeasible','Time-out')

    figure(3)
    scatter3(greenX',greenY',greenN',3,'g','.')
    hold all
    scatter3(redX',redY',redN',3,'r','.')
    scatter3(blackX',blackY',blackN',3,'black','.')
    xlabel('Tank Radius [m]')
    ylabel('Number of Tanks')
    zlabel('Optimal Num. Tanks')
    title('Starting Point Sensitivity Analysis, Number of Tanks')
    legend('Optimal','Infeasible','Time-out')
%}