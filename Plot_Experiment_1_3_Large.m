%Authors: Feng-Yi Liao & Yang Zheng
%         SOC Lab @UC San Diego
clc; clear; close;

%parameter setting 
width  = 8;     % Width in inches
height = 6;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 11;      % Fontsize
lw     = 1.5;      % LineWidth
msz    = 8;       % MarkerSize

filename = "n1000m200dr997";
filename_result = "results_rdSDPs\"+ filename + "_result.mat";
load(filename_result);
filename_solution = "examples\randomSDPs\" + filename+ ".mat";
load(filename_solution);

%range 
maxiter  = 300; 
xrange   = 1:maxiter;
xrange_p = xrange+1;

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



semilogy(abs((Out_Primal_0_2.Obj(xrange_p)-Optimal.Cost)/Optimal.Cost),"-.",'Color',COLORS(1,:),'LineWidth',lw);
hold on 
semilogy(abs((Out_Primal_0_3.Obj(xrange_p)-Optimal.Cost)/Optimal.Cost),"-.",'Color',COLORS(2,:),'LineWidth',lw); 
semilogy(abs((Out_Primal_0_4.Obj(xrange_p)-Optimal.Cost)/Optimal.Cost),"-.",'Color',COLORS(3,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_1_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'Color',COLORS(4,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_2_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'Color',COLORS(5,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_3_1.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"--",'Color',COLORS(6,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_0_2.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'Color',COLORS(7,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_0_3.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'Color',COLORS(8,:),'LineWidth',lw); 
semilogy(abs((-Out_Dual_0_4.Obj(xrange)-Optimal.Cost)/Optimal.Cost),"-",'Color',COLORS(9,:),'LineWidth',lw); 

legend('D(0,2)','D(0,3)','D(0,4)','P(1,1)','P(2,1)','P(3,1)',...
       'P(0,2)','P(0,3)','P(0,4)','NumColumns',3);
xlabel('Iteration','interpreter','latex');
%ylabel('$\epsilon_{\mathrm{o}}$','interpreter','latex');
ylabel('Cost Opt.','interpreter','latex');
title('rank$(X^\star)=3$','interpreter','latex');
set(gca,'FontSize',14);
ax = gca;
ax.YLabel.FontSize = 24;
ylim([5*10^-8,10^6]);
print("results_rdSDPs\LowRank-1",'-depsc','-tiff');

normb = norm(b_sdp);
normc = norm(c_sdp);
fprintf('name | e_p | e_d | e_r | e_g | e_o\n');
fprintf('Out_Primal_0_2   %7.2e   %7.2e  %7.2e   %7.2f\n',Out_Primal_0_2.DescentPrimalSemiFeasi(end),Out_Primal_0_2.DescentRelativeDFeasi(end),Out_Primal_0_2.DescentRelativeGap(end),abs((Out_Primal_0_2.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Primal_0_3   %7.2e   %7.2e  %7.2e   %7.2f\n',Out_Primal_0_3.DescentPrimalSemiFeasi(end),Out_Primal_0_3.DescentRelativeDFeasi(end),Out_Primal_0_3.DescentRelativeGap(end),abs((Out_Primal_0_3.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Primal_0_4   %7.2e   %7.2e  %7.2e   %7.2f\n',Out_Primal_0_4.DescentPrimalSemiFeasi(end),Out_Primal_0_4.DescentRelativeDFeasi(end),Out_Primal_0_4.DescentRelativeGap(end),abs((Out_Primal_0_4.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_1_1  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_1_1.DescentDualSemiFeasi(end),Out_Dual_1_1.DescentRelativePFeasi(end),Out_Dual_1_1.DescentRelativeGap(end),abs((-Out_Dual_1_1.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_2_1  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_2_1.DescentDualSemiFeasi(end),Out_Dual_2_1.DescentRelativePFeasi(end),Out_Dual_2_1.DescentRelativeGap(end),abs((-Out_Dual_2_1.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_3_1  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_3_1.DescentDualSemiFeasi(end),Out_Dual_3_1.DescentRelativePFeasi(end),Out_Dual_3_1.DescentRelativeGap(end),abs((-Out_Dual_3_1.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_0_2  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_0_2.DescentDualSemiFeasi(end),Out_Dual_0_2.DescentRelativePFeasi(end),Out_Dual_0_2.DescentRelativeGap(end),abs((-Out_Dual_0_2.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_0_3  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_0_3.DescentDualSemiFeasi(end),Out_Dual_0_3.DescentRelativePFeasi(end),Out_Dual_0_3.DescentRelativeGap(end),abs((-Out_Dual_0_3.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
fprintf('Out_Dual_0_4  %7.2e   %7.2e  %7.2e   %7.2e\n',Out_Dual_0_4.DescentDualSemiFeasi(end),Out_Dual_0_4.DescentRelativePFeasi(end),Out_Dual_0_4.DescentRelativeGap(end),abs((-Out_Dual_0_4.DescentCost(end)-Optimal.Cost)/Optimal.Cost));
