for i=1:1001
    ar(i) = i; 
    rm(i)= sqrt(1/ar(i));   
end 
arr_111=[32 20 16 10 8 5 4 2 1];
arr_211=[32 20 16 10 8];
arr_333=[8];
rmc_111=[0.1753 0.2354 0.2695 0.3385 0.4275 0.7805 0.9785 4.425 56.52];
rmc_211=[0.195  0.245   0.29    0.358   4.275000e-01];
rmc_333=[1.2];
fontt=22;

subplot(2,1,1);
semilogy(ar, rm,'k:','LineWidth',2);
title('$\textit{Critical }R_{mc}$','interpreter','Latex', 'fontsize',20)
ylabel('$R_{mc}$','interpreter','Latex', 'fontsize',fontt);
hold on;
plot(arr_111,rmc_111, 'ro','LineWidth',1.5);
plot(arr_211,rmc_211, 'rx','LineWidth',1.5);
plot(arr_333,rmc_333, 'bx','LineWidth',1.5);
hold off;
xlim([5 35]);
h_legend=legend('$\textit{Asymptotic Approximation}$', '${\bf{u}}_{ABC}, ar=L\times1\times1$','${\bf{u}}_{ABC}, ar=L\times2\times1$','${\bf{u}}_{III}, ar=L\times2\times1$');
set(h_legend,'interpreter','Latex','FontSize',14,'box', 'off')

subplot(2,1,2);
for ii=1:size(arr_111,2)
    for jj=1:size(ar,2)
        if ar(jj)==arr_111(ii)
            error_111(ii)= abs(rm(jj)-rmc_111(ii))/rm(jj);
        end 
    end
end
for ii=1:size(arr_211,2)
    for jj=1:size(ar,2)
        if ar(jj)==arr_211(ii)
            error_211(ii)= abs(rm(jj)-rmc_211(ii))/rm(jj);
        end 
    end
end

semilogy(arr_111,error_111, 'ko','LineWidth',1.5);
hold on;
semilogy(arr_211,error_211, 'kx','LineWidth',1.5);
hold off;
xlim([5 35]);
ylabel('$\% Error$','interpreter','Latex', 'fontsize',fontt);
xlabel('$\textit{Box Aspect Ratio}$','interpreter','Latex', 'fontsize',fontt);

% subplot(3,1,3);
% plot(ar, rm,'k-','LineWidth',2);
% xlabel('$\textit{Aspect Ratio}$','interpreter','Latex', 'fontsize',fontt);
% ylabel('$R_{mc}$','interpreter','Latex', 'fontsize',fontt);
% ylim([0.03 0.4]);
% xlim([1 50]);
% str = '$$ R_{mc}=\left(\frac{1}{AR}\right)^{1/2} $$';
% text(33.3,0.23,str,'Interpreter','latex','fontsize',16)
% hold on;
% 
% for jj=1:6
%     colorr = rand(1,3);
%     plot([1 arr_111(jj)],[sqrt(1/arr_111(jj)) sqrt(1/arr_111(jj))], 'k:','LineWidth',1.5,'color',colorr);
%     plot([arr_111(jj) arr_111(jj)],[0 sqrt(1/arr_111(jj))], 'k:','LineWidth',1.5, 'color',colorr);
% end 
% hold off;