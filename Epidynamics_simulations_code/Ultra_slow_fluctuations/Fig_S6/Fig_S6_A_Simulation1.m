%% Status epilepticus below SH

% % create the vector Status_def from the 4th (i=4) simulation of RepRep2=4

% t=size(timeseriesrep,2);
% v=size(timeseriesrep,3);
% X=zeros(t,v);
% 
% for i=1:t
%     for j=1:v
%         X(i,j)=timeseriesrep(4,i,j);
%     end
% end

%%

load('Status_beforeSH.mat')
%%

% plot curve followed in parameter space
openfig('../../0_COD3MODEL/FOCUS_r0p4/Unfolding_r0p4.fig')
hold on

    Au=R*(G'*cos(X(:,4))'+H'*sin(X(:,4))'); % A(u) and B(u)
    Bu=R*(L'*cos(X(:,5))'+M'*sin(X(:,5))');
    
    Eu=Au/R;                        % vectors for path from A(u) to B(u)
    Fu=cross(cross(Au,Bu),Au);
    for i=1:size(X(:,4),1)
    Fu(:,i)=Fu(:,i)/norm(Fu(:,i));
    end
   
    mu2=R*(Eu(1,:).*cos(X(:,3))'+Fu(1,:).*sin(X(:,3))'); 
    mu1=-R*(Eu(2,:).*cos(X(:,3))'+Fu(2,:).*sin(X(:,3))');
    nu=R*(Eu(3,:).*cos(X(:,3))'+Fu(3,:).*sin(X(:,3))');

dt=0.002;
tmin_plot=3000/dt;
tmax_plot=6000/dt;

plot3(mu2(tmin_plot:tmax_plot),-mu1(tmin_plot:tmax_plot),nu(tmin_plot:tmax_plot),'k','LineWidth',1.5)
%plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
set(gcf, 'Color', 'w')
hold off
%% plot flat

ph = atan2(-mu1,mu2);
th = acos(nu./(sqrt(mu2.^2+mu1.^2+nu.^2)));

% plot freq

openfig('../../0_COD3MODEL/Figures/Freq_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)
    %title('Frequency')
    xlabel('\phi')
    ylabel('\theta')
    
    % plot ampl

openfig('../../0_COD3MODEL/Figures/Ampl_Equazim.fig')
hold on
topolegend=[1,90,0];
axesm('MapProjection','eqaazim','Grid','on')
   
    path_lat=(pi/2-th)*180/pi;
    path_long=ph*180/pi;
    geoshow(path_lat(tmin_plot:tmax_plot),path_long(tmin_plot:tmax_plot),'Color','w','LineWidth',2)
    plotm(path_lat(tmin_plot),path_long(tmin_plot),'^r','LineWidth',0.5,'MarkerFaceColor','r','MarkerSize',6)
    plotm(path_lat(tmax_plot),path_long(tmax_plot),'xr','LineWidth',2,'MarkerSize',8)

    %title('Amplitude')
    xlabel('\phi')
    ylabel('\theta')