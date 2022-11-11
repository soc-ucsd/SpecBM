clc; clear; close;
filename_result = 'examples\Result\n100m100dr97\n100m100dr97-result-test5.mat';
load(filename_result);
filename_solution = 'examples\n100m100dr97';
load(filename_solution);

%parameter setting 
width  = 8;     % Width in inches
height = 6;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 11;      % Fontsize
lw     = 1.5;      % LineWidth
msz    = 8;       % MarkerSize

%range 
xrange = 1:200;

%for plot
res         = 10;
xrange_plot = 1:res:200;
N           = 9;%number of colors
C           = linspecer(N);
     
COLORS = [0.4660, 0.6740, 0.1880;
    0.3010, 0.7450, 0.9330;
    0.6350, 0.0780, 0.1840;
    0, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.4940, 0.1840, 0.5560;
    1, 0, 0;
    0, 1, 0;
    0, 0, 1];

figure();
pos = get(gcf, 'Position');
%set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

semilogy(abs((Out_Primal_0_3.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-.",'LineWidth',lw,'Color',COLORS(1,:)); 
hold on 
semilogy(abs((Out_Primal_0_5.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-.",'LineWidth',lw,'Color',COLORS(2,:)); 
semilogy(abs((Out_Primal_0_10.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-.",'LineWidth',lw,'Color',COLORS(3,:)); 

semilogy(abs((-Out_Dual_2_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'LineWidth',lw,'Color',COLORS(4,:));
semilogy(abs((-Out_Dual_4_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'LineWidth',lw,'Color',COLORS(5,:)); 
semilogy(abs((-Out_Dual_9_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'LineWidth',lw,'Color',COLORS(6,:)); 

semilogy(abs((-Out_Dual_0_3.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'LineWidth',lw,'Color',COLORS(7,:)); 
semilogy(abs((-Out_Dual_0_5.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'LineWidth',lw,'Color',COLORS(8,:)); 
semilogy(abs((-Out_Dual_0_10.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'LineWidth',lw,'Color',COLORS(9,:)); 

legend('P(0,3)','P(0,5)','P(0,10)','D(2,1)','D(4,1)','D(9,1)',...
       'D(0,3)','D(0,5)','D(0,10)','NumColumns',3);
xlabel('Iteration','interpreter','latex');
ylabel('$\frac{F-F^\star}{{|F^\star|}}$','interpreter','latex');
title('$X_\star$ is low rank','interpreter','latex');
set(gca,'FontSize',14);
ax = gca;
ax.YLabel.FontSize = 24;
ylim([10^-9,10^4]);
print("LowRank",'-depsc','-tiff');