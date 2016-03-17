
subplot(3,1,1);
plot( box32(:,1), box32(:,2),'k:',     box32(:,1),box32(:,16),'-' ,      box32(:,1),box32(:,17),'-' ,      box32(:,1),box32(:,18),'-' ,      box32(:,1),box32(:,19),'-' ,      box32(:,1),box32(:,20),'-' ,      box32(:,1),box32(:,21),'-' ,      box32(:,1),box32(:,22),'-' ,      box32(:,1),box32(:,23),'-' ,      box32(:,1),box32(:,24),'-' ,      box32(:,1),box32(:,25),'-' ,      box32(:,1),box32(:,26),'-',      box32(:,1),box32(:,27),'-',box32(:,1),box32(:,28),'-');
title({'Box of AR=32 for the modified ABC flow', 'Growth rate, p(R_m) from the dispersion relation'; },'FontSize',15);
xlabel('x/L_x')
leg1='$$\cos^2(\frac{2\pi}{L_x}x)$$';
leg2='$$R_m=0.003$$';
leg3='$$R_m=0.004$$';	
leg4='$$R_m=0.005$$';
leg5='$$R_m=0.01$$';
leg6='$$R_m=0.03$$';	
leg7='$$R_m=0.035$$'; 
leg8='$$R_m=0.04$$'; 
leg9='$$R_m=0.045$$';	
leg10='$$R_m=0.05$$';
leg11='$$R_m=0.055$$';
leg12='$$R_m=0.06$$'; 
leg13='$$R_m=0.065$$';
leg14='$$R_m=0.1$$';

h =legend(leg1,leg2, leg3,leg4,leg5,leg6,leg7,leg8,leg9,leg10,leg11,leg12,leg13, leg14);
set(h,'Interpreter','latex');
set(h,'FontSize',15);

subplot(3,1,2);
plot( box32(:,1), box32(:,2),'k:',box32(:,1),box32(:,3), '-',       box32(:,1),box32(:,3),'-' ,      box32(:,1),box32(:,4),'-' ,      box32(:,1),box32(:,5),'-' ,      box32(:,1),box32(:,6),'-' ,      box32(:,1),box32(:,7),'-' ,      box32(:,1),box32(:,8),'-' ,      box32(:,1),box32(:,9),'-' ,      box32(:,1),box32(:,10),'-' ,      box32(:,1),box32(:,11),'-' ,      box32(:,1),box32(:,12),'-' ,      box32(:,1),box32(:,13),'-' ,      box32(:,1),box32(:,14),'-');
xlabel('x/L_x')
title('% Difference of growth rates, p(R_m=0.1) - p(R_m)','FontSize',15);
ylim([0 100])


subplot(3,1,3);
plot( box32(:,1), box32(:,2),'k:',box32(:,1),box32(:,3), '-',       box32(:,1),box32(:,3),'-' ,      box32(:,1),box32(:,4),'-' ,      box32(:,1),box32(:,5),'-' ,      box32(:,1),box32(:,6),'-' ,      box32(:,1),box32(:,7),'-' ,      box32(:,1),box32(:,8),'-' ,      box32(:,1),box32(:,9),'-' ,      box32(:,1),box32(:,10),'-' ,      box32(:,1),box32(:,11),'-' ,      box32(:,1),box32(:,12),'-' ,      box32(:,1),box32(:,13),'-' ,      box32(:,1),box32(:,14),'-');
xlabel('x/L_x')
ylim([0 2])

figure(2)
subplot(2,1,1);
plot( rates(:,1), rates(:,2),'k:',     rates(:,1),rates(:,3),'-' ,rates(:,1),rates(:,4),'-' ,rates(:,1),rates(:,5),'-' ,rates(:,1),rates(:,6),'-' );
leg1='$$\cos^2(\frac{2\pi}{L_x}x)$$';
leg2='$$AR=32$$';
leg3='$$AR=20$$';	
leg4='$$AR=16$$';
leg5='$$AR=10$$';
h =legend(leg1,leg2, leg3,leg4,leg5,'Location','eastoutside');
set(h,'Interpreter','latex');
set(h,'FontSize',15);
title({'Growth rate, p(R_m) from the dispersion relation','For the modified ABC flow at R_m=0.1'  },'FontSize',15);
xlabel('x/L_x')

subplot(2,1,2);
plot( rates(:,1), rates(:,2),'k:',     rates(:,1),rates(:,7),'-' ,rates(:,1),rates(:,8),'-' ,rates(:,1),rates(:,9),'-' ,rates(:,1),rates(:,10),'-' );
leg1='$$\cos^2(\frac{2\pi}{L_x}x)$$';
leg2='$$AR=32$$';
leg3='$$AR=20$$';	
leg4='$$AR=16$$';
leg5='$$AR=10$$';
h =legend(leg1,leg2, leg3,leg4,leg5,'Location','eastoutside' );
set(h,'Interpreter','latex');
set(h,'FontSize',15);
title({'For the modified ABC flow at R_m=0.7' },'FontSize',15);
xlabel('x/L_x')