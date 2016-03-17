%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PLOT %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boxes=[32,20,16,10,8,5,4,2,1];
   

figure(1)
subplot(2,1,1)
hold on;
for i=2:size(boxes,2)
    % Generate random colour
    colorr = rand(1,3);
    
    % Eliminate zero entries 
    for k=1:size(rm_vs_box_data,1)
        kk=1;
        if(rm_vs_box_data(k,i) ~= 0.0)
            rm(kk,i)= rm_vs_box_data(k,1);
            rates(kk,i)= rm_vs_box_data(k,i);
            kk=kk+1;
        end    
    end
    plot( rm(:,i),rates(:,i),'--*', 'color', colorr);
   
end 
hold off;

title({'R_m vs Growth Rate for boxes of aspect ratio X:1:1','for 111-ABC flow with k^u=1'});
ylabel('Growth Rate');
xlim([0,130]);
ylim([0,0.2]);

subplot(2,1,2)
hold on;
for i=2:size(boxes,2)
    colorr = rand(1,3);
    plot( rm_vs_box_data(:,1),rm_vs_box_data(:,i),'--*', 'color', colorr);
end 
hold off;
ylabel('Growth Rate');
xlabel('R_m');
xlim([20,50]);