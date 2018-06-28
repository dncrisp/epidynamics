% Plots the unfolding of the deg. TB singularity on a sphere. See Saggio et
% al. 2017 JMN for details

% load bifurcation curves obtained with Matcont and CL_Matcont
load('WS_FOCUS_r0p4_bifcurves_and_statespace.mat')

figure

% curves appearence settings
foldcolor=[1, 150/255, 20/255];  % orange
hopfcolor=[0, 204/255, 0];       % green
homcolor=[0, 153/255, 1];        % azzurro
fccolor=[1, 0,153/255];    % magenta
bifwidth=1.5;
bifcolor=[foldcolor; hopfcolor; homcolor; fccolor]


for index_plot = 1
    try
        % plot bif sets in parameter space
        
        hold on
        nmax=120;
        p3=plot3(HomTB1{index_plot,2}(3,1:nmax), HomTB1{index_plot,2}(4,1:nmax), HomTB1{index_plot,2}(5,1:nmax),'Color',bifcolor(3,:),'LineWidth',bifwidth)
        nmax=78;
        plot3(HomTB2{index_plot,2}(3,1:nmax), HomTB2{index_plot,2}(4,1:nmax), HomTB2{index_plot,2}(5,1:nmax),'Color',bifcolor(3,:),'LineWidth',bifwidth)
        nmax=130;
        plot3(HomMiddle{index_plot,2}(3,1:nmax), HomMiddle{index_plot,2}(4,1:nmax), HomMiddle{index_plot,2}(5,1:nmax),'Color',bifcolor(3,:),'LineWidth',bifwidth)
        nmax=130;
        plot3(HomMiddleb{index_plot,2}(3,1:nmax), HomMiddleb{index_plot,2}(4,1:nmax), HomMiddleb{index_plot,2}(5,1:nmax),'Color',bifcolor(3,:),'LineWidth',bifwidth)
        
        p1=plot3(fold{index_plot,2}(3,:), fold{index_plot,2}(4,:), fold{index_plot,2}(5,:),'Color',bifcolor(1,:),'LineWidth',bifwidth)
        plot3(hopf_curve{index_plot,2}(3,:), hopf_curve{index_plot,2}(4,:), hopf_curve{index_plot,2}(5,:),'Color',bifcolor(2,:),'LineWidth',bifwidth)
        
        p2=plot3(hopf_curve{index_plot,2}(3,:), hopf_curve{index_plot,2}(4,:), hopf_curve{index_plot,2}(5,:),'Color',bifcolor(2,:),'LineWidth',bifwidth)
        
        for i=1:3
            p4=plot3(fold_cycle{index_plot,2,i}(3,:), fold_cycle{index_plot,2,i}(4,:), fold_cycle{index_plot,2,i}(5,:),'Color',bifcolor(4,:),'LineWidth',bifwidth)
            % plot3(fold_cycle_b{index_plot,2,i}(3,:), fold_cycle_b{index_plot,2,i}(4,:), fold_cycle_b{index_plot,2,i}(5,:),'b')
        end
        
        % title('Cod 3 Focus type')
        xlabel('\mu_2')
        ylabel('-\mu_1')
        zlabel('\nu')
        legend([p1,p2,p3,p4],'Saddle-Node','Hopf','Saddle-Homoclinic','Fold Limit Cycle')
    end
end
hold on
[x y z]=sphere(60);
r=0.399;
surf(r*x,r*y,r*z)
colormap([206 224 248]/255)
shading interp;
alpha(0.6)
axis equal
%title('Unfolding')
set(gcf,'color','w');