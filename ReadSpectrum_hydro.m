close all
% Spectrum of flow u in a box [10,2,1] from the hydro sims
% 
box=10;
% What files are out there
rep=['/data/novadisk/vs391/hydro/u_iii/'];
%rep=['/data/novadisk/vs391/hydro/u_abc/box_10'];
infofile= [rep,'/info'];
% Locate these files
rep1=['/data/novadisk/vs391/hydro/u_iii/run_'];
%rep1=['/data/novadisk/vs391/hydro/u_abc/box_10/run_'];

% Read from info file what Re's we have
% and the corresponding location
full_file=importdata(infofile);
timevar=full_file;%.data;
tblA = table(timevar(:,1),timevar(:,2), timevar(:,3));
% Sort the rows of the table based on Rm
tblB = sortrows(tblA,3); 
run = tblB{1:end,1}; 
cas = tblB{1:end,2}; 
Re  = tblB{1:end,3}; 

% Assign to the files to read out in the loop
for j=1:size(Re)
    files1{j}=[rep1,num2str(run(j)),'/spectrum',num2str(cas(j)),'.dat'];
end
    color{2} = [1,0,0];
    color{1} = [0,0,0];
mean_start=1;
for ii=1:size(Re)
    %color{ii}=rand(1,3); 
    
    nspec=21;    % 9 for old spectrum (no transfer). 15 for recent versions, 21 for v5.0 version
    spectruml=importdata(files1{1,ii});
    spectrum.k=spectruml(1,1:(end-1));
    spectrum.n=spectruml(2,1:(end-1));
    spectrum.vx=spectruml(3:nspec:end,2:end);
    spectrum.vy=spectruml(4:nspec:end,2:end);
    spectrum.vz=spectruml(5:nspec:end,2:end);
    spectrum.bx=spectruml(6:nspec:end,2:end);
    spectrum.by=spectruml(7:nspec:end,2:end);
    spectrum.bz=spectruml(8:nspec:end,2:end);
    spectrum.th=spectruml(9:nspec:end,2:end);
    spectrum.vxvy=spectruml(10:nspec:end,2:end);
    spectrum.bxby=spectruml(11:nspec:end,2:end);
    spectrum.ad_vx=spectruml(12:nspec:end,2:end);
    spectrum.ad_vy=spectruml(13:nspec:end,2:end);
    spectrum.ad_vz=spectruml(14:nspec:end,2:end);
    spectrum.tr_bx=spectruml(15:nspec:end,2:end);
    spectrum.tr_by=spectruml(16:nspec:end,2:end);
    spectrum.tr_bz=spectruml(17:nspec:end,2:end);
    spectrum.tr_vx=spectruml(18:nspec:end,2:end);
    spectrum.tr_vy=spectruml(19:nspec:end,2:end);
    spectrum.tr_vz=spectruml(20:nspec:end,2:end);
    spectrum.hel=spectruml(21:nspec:end,2:end)+spectruml(22:nspec:end,2:end)+spectruml(23:nspec:end,2:end);


    spectrum.t=transpose(spectruml(4:nspec:end,1));

    clear spectrum1;

    % Compute the fluxes
    flux.adk=-cumsum(mean(spectrum.ad_vx(mean_start:end,:)+spectrum.ad_vy(mean_start:end,:)+spectrum.ad_vz(mean_start:end,:)));
    flux.adm=-cumsum(mean(spectrum.tr_vx(mean_start:end,:)+spectrum.tr_vy(mean_start:end,:)+spectrum.tr_vz(mean_start:end,:)+spectrum.tr_bx(mean_start:end,:)+spectrum.tr_by(mean_start:end,:)+spectrum.tr_bz(mean_start:end,:)));

    %%%%%%%%%
    % PLOTS %
    %%%%%%%%%
        hFig = figure(box);
        set(hFig, 'Position', [100, 60, 1049, 900]);
        
        subplot(3,4,ii);  
        set(gca, 'FontSize', 12)    
            h12=loglog(spectrum.k,(spectrum.vx(1,:)+spectrum.by(1,:)+spectrum.vz(1,:))/max(spectrum.vx(1,:)+spectrum.vy(1,:)+spectrum.vz(1,:)),'o',...
                            'LineStyle', 'none',...
                            'color',color{2},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{2}, ...
                            'MarkerFaceColor', color{2}, ...
                            'MarkerSize',4.5);       
                        hold on;
                        
            h13=loglog(spectrum.k,(spectrum.vx(end,:)+spectrum.vy(end,:)+spectrum.vz(end,:))/max(spectrum.vx(end,:)+spectrum.vy(end,:)+spectrum.vz(end,:)),'o',...
                            'LineStyle', 'none',...
                            'color',color{1},...
                            'LineWidth',1.5,...
                            'MarkerEdgeColor', color{1}, ...
                            'MarkerFaceColor', color{1}, ...
                            'MarkerSize',4.5);              
            title(['R_e=',num2str(Re(ii))]);
            ylim([1e-3 1])
            xlim([1e-2 1e1])
            if ii > size(Re) -4 
                xlabel('k','fontsize',16, 'Interpreter', 'latex');
            end
            
            plot([1/box 1/box],[1e-10 1],'k:');
            plot([1-1/box 1-1/box],[1e-10 1],'k:');
            plot([1 1],[1e-10 1],'k:');
            plot([1+1/box 1+1/box],[1e-10 1],'k:');
            if ii==size(Re)-3
                legend([h12 h13],{'initial','final'},'Location','southwest');
            end
end
