

clear all;

%------- SETTINGS - INTEGRATION
options = odeset('MaxStep', 0.1);   % maximal time integration step

fs_exp = 512;                 %[fs]=1/s %system sampling frequency
dt_exp = 1 / fs_exp;               %[dt]=s

dt=dt_exp;
tmax= 20000;
tspan = (0 : dt : tmax)';    % time span to integrate over
xinit=[0;0;0;0;0];            % initial conditions (must be a column)

%------- SETTINGS - MODEL
b = 1.0;                            % focus
R = 0.4;                            % radius in unfolding
dstar =0.3 ;                        % threshold for slow variable inversion
N_equil=1;                          % resting state is upper branch of equilibria

%noise
dim = [1,2,3,4,5]; % additive noise to variable in dim

% settings used to have status epilepticus;
c = 0.0001; % velocity of slow variable
cA = 0; % velocity of modulation of A
cB=0; % velocity of modulation of B
Ain=[0.3483,0.03698,0.1931];
Bin=[0.3331,0.074,0.2087];
Aend=[0.279,0.2187,0.1854];
Bend=Aend;

%------- INTEGRATION

% vectors G,H,M,L to describe path
G=Ain/R;
H=cross(cross(Ain,Aend),Ain);
H=H/norm(H);

L=Bin/R;
M=cross(cross(Bin,Bend),Bin);
M=M/norm(M);

% noise scaling factor
a=[0.0005, 0.0002, 0.0001, 0.00005];

for reprep=1:length(a)
    
    reprep
    % noise variance for each variable
    sigma=a(reprep)*[10,10,1,5,5];
    % repetion of simulation
    Nrep=5;
    Nt=length(tspan);
    timeseriesrep=zeros(Nrep,Nt,5);
    
    for k=1:Nrep
        
        k
        clear Fx0 x1 Fx1 nx
        %-------------------------------------------------------
        
        
        N = 1:nnz(tspan);
        X = zeros(5,length(tspan));
        x = xinit;
        r = zeros(1,length(N));

        noisevar=length(dim);
        for i=1:noisevar
            Rn(dim(i),:) = sigma(i)*randn(1,length(N));
        end
        
        for n = N
            r = Rn(:,n);
            
            % Euler-Meruyama method
            Fx0 = Transitions_Cod3_Model(tspan,x,options,b,c,R,G,H,L,M,dstar,cA,cB,N_equil);
            nx = x + dt*Fx0 + sqrt(dt)*r;
            x=nx;
            X(:,n) = x;
        end
        
        X=X';
        
        timeseriesrep(k,:,:)=X;
    end      

title=sprintf('RepRep2=%d', reprep);
save(title)
    end
    
    %% plot curve followed in parameter space
    for i=1:Nrep
        
        X = zeros(length(tspan),5);
        for ix=1:5
            for jx=1:length(tspan)
                X(jx,ix)=timeseriesrep(i,jx,ix);
            end
        end
        
        openfig('../../../0_COD3MODEL/FOCUS_r0p4/Unfolding_r0p4.fig')
        hold on
        
        Au=R*(G'*cos(X(:,4))'+H'*sin(X(:,4))'); % A(u) and B(u)
        Bu=R*(L'*cos(X(:,5))'+M'*sin(X(:,5))');
        
        Eu=Au/R;                        % vectors for path from A(u) to B(u)
        Fu=cross(cross(Au,Bu),Au);
        for ii=1:size(X(:,4),1)
            Fu(:,ii)=Fu(:,ii)/norm(Fu(:,ii));
        end
        
        mu2=R*(Eu(1,:).*cos(X(:,3))'+Fu(1,:).*sin(X(:,3))');
        mu1=-R*(Eu(2,:).*cos(X(:,3))'+Fu(2,:).*sin(X(:,3))');
        nu=R*(Eu(3,:).*cos(X(:,3))'+Fu(3,:).*sin(X(:,3))');
        
        %tmin_plot=5000/dt; tmax_plot=6000/dt;
        %plot3(mu2(tmin_plot:tmax_plot),-mu1(tmin_plot:tmax_plot),nu(tmin_plot:tmax_plot),'r','LineWidth',1.5)
        plot3(mu2,-mu1,nu,'r','LineWidth',1.5)
        set(gcf, 'Color', 'w')
        %title(str)
        hold off
   
   
    
    
    figure
    newblue=[0, 113/255, 188/255];
    newred=[216/255, 82/255, 24/255];
    plot(tspan, X(:,1),'Color',newblue);
    hold on
    plot( tspan,X(:,3),'Color',newred);
    plot( tspan,X(:,4));
    % hold on
    % plot(tspan, Rn,'Color','g');
    hold off
    %title(str)
    xlabel('time')
    ylabel('x')
    legend('x','z')
    set(gcf, 'Color', 'w')
    
    end



