%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	This file is part of the Snoopy code.

%   Snoopy code is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    Snoopy code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Read spectrum written by latest versions of Snoopy (this is a beta
% version).
files1=[ {'/data/novadisk/vs391/hydro/hydro_lx10_re12/spectrum.dat'}...  
         {'/data/novadisk/vs391/snoopy_kinematic_dynamo/u_iii/snoopy1/spectrum.dat'}];
legendInfo1 = [     {'u_{abc} in [10,2,1]'}...
                    {'u_{III} in [10,2,1]'}];
% Are we plotting B and U - set to 1 if yes? 
B_flag = 0;
u_flag = 1;
%rep='/Volumes/Idris_Workdir/snoopy/mri/768x384x192/Pm1/';

mean_start=1;
for ii=1:size(files1,2)
    %color{ii}=rand(1,3); 
    color{1} = [1,0,0];
    color{2} = [0,0,0];
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


    %figure(1)
    %loglog(spectrum.k,spectrum.n);
    %title('Number of modes');
    %hold on;
    
    figure(2)
    if (B_flag==1 && u_flag==1)
        loglog(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)),spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)),spectrum.k,1e-1*spectrum.k.^(-5/3));
        title('Spectrum');
        legend('Kinetic','Magnetic','K43')
    elseif (B_flag==1 && u_flag==0)
        loglog(spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)));
        title('Magnetic Spectrum');
        legend(legendInfo1,'Location','east')
    elseif (B_flag==0 && u_flag==1)
        semilogy(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)) ...
            ,':o','LineWidth',1.5, ...
            'color', color{ii}, ...
            'MarkerEdgeColor', color{ii}, ...
            'MarkerFaceColor', color{ii}, ...
            'MarkerSize',3.5);
        title('Kinetic Energy Spectrum');
        xlabel('$k [2\pi L^{-1}]$','fontsize',16, 'Interpreter', 'latex');
        ylabel('$E_K [L^2T_{turnover}^{-2}]$','fontsize',16, 'Interpreter', 'latex'  );
        legend(legendInfo1,'Location','northeast')
    end
    hold on;
if (1==2)
    figure(3)
    if (B_flag==1 && u_flag==1)
        semilogx(spectrum.k,mean(spectrum.vxvy(mean_start:end,:)),spectrum.k,mean(spectrum.bxby(mean_start:end,:)));
        title('Transport');
        legend('Reynolds','Maxwell')
    elseif(B_flag==1 && u_flag==0)
        semilogx(spectrum.k,mean(spectrum.bxby(mean_start:end,:)));
        title('Transport');
        legend('Maxwell')
    elseif(B_flag==0 && u_flag==1)
        semilogx(spectrum.k,mean(spectrum.vxvy(mean_start:end,:)));
        title('Transport');
        legend('Reynolds')
    end
    hold on;
    
    if (u_flag==1)
    figure(4)
        loglog(spectrum.k,spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:))
        title('Kinetic spectrum in time')
    end


    if (B_flag==1)
        figure(5)
        loglog(spectrum.k,spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))
        title('Magnetic spectrum in time')
    end
    hold on;
    figure(6)
    if (B_flag==1 && u_flag==1)
        semilogx(spectrum.k,flux.adk,spectrum.k,flux.adm);
        title('Fluxes');
        legend('Kinetic transfer','Magnetic transfer')
    elseif(B_flag==1 && u_flag==0)
        semilogx(spectrum.k,flux.adm);
        title('Fluxes');
        legend('Magnetic transfer')
    elseif(B_flag==0 && u_flag==1)
        semilogx(spectrum.k,flux.adk);
        title('Fluxes');
        legend('Kinetic transfer')
    end
    hold on;

    figure(7)
    if (B_flag==1 && u_flag==1)
        loglog(spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)),spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)),spectrum.k,1e-1*spectrum.k.^(1/3));
        title('dissipation');
        legend('Kinetic','Magnetic','K43')
    elseif(B_flag==1 && u_flag==0)
        loglog(spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)));
        title('dissipation');
        legend('Magnetic')
    elseif(B_flag==0 && u_flag==1)
        semilogy(spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)));
        title('dissipation');
        legend('Kinetic')
    end
    hold on;

    figure(8)
    semilogx(spectrum.k,spectrum.k.*mean(spectrum.hel(mean_start:end,:))),
    title('Kinetic helicity');
    hold on;

    figure(9)
    semilogx(spectrum.k,spectrum.hel(mean_start:end,:));
    title('Kinetic helicity in time')
    hold on;
end
end