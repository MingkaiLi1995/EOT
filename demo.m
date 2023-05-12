%*********************************************************************
%**    A real data implementation of the extended object tracking   **
%**    algorithms based on the paper                                **
%**    "Tracking of Elliptical Object with Unknown                  **
%**     but Fixed Lengths of Axes"                                  **
%**    IEEE Transactions on Aerospace and Electronic Systems.       **
%**    Mingkai Li, Jian Lan, and X. Rong Li                         **
%**                                                                 **
%*********************************************************************

clc
clear
close all
%% load data

set(0,'defaulttextinterpreter','tex')
addpath sample_function
addpath realdata

load 'point_measurements.dat' % point measurements
load 'vertices_of_groundtruth.dat' % groundtruth: a 12.97m*2.74m bus
load 'number_of_point_measurements_per_scan.dat' % number of point measurements
%% Paramaters

time_steps = 70;
T = 0.5;

% object parameters

scale1.mean = 2/3;
scale1.variance = 1/18;
%% Process Model

n_orientation = 2;
F_k_cv = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
F_k_extent = eye(n_orientation);
F_k = blkdiag(F_k_extent,F_k_cv);

C_w_p = blkdiag(0.1^2,0.1^2);   
sigmaQ = sqrt(5);
C_w_r = sigmaQ^2 * [T^3/3 T^2/2 0 0
                    T^2/2   T   0 0 
                    0 0 T^3/3 T^2/2
                    0 0 T^2/2 T];
C_w = blkdiag(C_w_p,C_w_r);
C_v = diag([0.25^2, 0.25^2]);

HH = [0 0 1 0 0 0
      0 0 0 0 1 0];
%% setting parameters prior

x_srt1 = zeros(n_orientation + 4, time_steps);
PP = zeros((n_orientation + 4),(n_orientation + 4),(time_steps));
Nsum = 0;
for t = 1:time_steps
	%% get measurements
 
    N = number_of_point_measurements_per_scan(t);
    disp([ 'time step:' num2str(t) ', ' num2str(N) ' measurements']);
    y = point_measurements(:, Nsum + 1: Nsum + N);
    Nsum = Nsum + N;
    Y(1,t) = {y};
	%% Visulize groundtruth and measurements
    
    if mod(t,7) == 0
        realX = vertices_of_groundtruth(3*t-2:3*t-2,:);
        realY = vertices_of_groundtruth(3*t-1:3*t-1,:);
        hold on
        gt_plot = plot([realX realX(1)],[realY realY(1)],'color','k');
        h_object = fill([realX realX(1)], [realY realY(1)], [0.98,0.9,0.1],'edgealpha',0);
        axis equal
        drawnow;
        y = Y{t};
        N = size(y,2);
        for j = 1 : N
            newMeasurement = y(:,j);
            drawnow;
            plot(newMeasurement(1), newMeasurement(2), 'k.','lineWidth', 0.5);
            hold on
        end
        hold on
    end
end
    
 %% EOT-EM
    tol = 1e-6;
    maxIter = 10;
    l=[5 1];
    Q = -inf(1,maxIter);
for iter = 2 : maxIter
    
    disp([num2str(iter-1) ' EOT-EM Step']);
    %% initialization for realdata
    
    x0 = [0.0248  0.9997 1062 -0.5 853 1]';
    P0 = blkdiag(0.04*diag(ones(1, n_orientation)),9,1,9,1);
    x_srt1(:,1) = x0;
    C_x1 = P0;
    %% state estimation
    
     for t = 1:time_steps
        y = Y{t};
        N = size(y,2);
        if l(2)>l(1)
           l_m = l(2);l(2) = l(1);l(1) = l_m;
        end
        for j = 1 : N
            newMeasurement = y(:,j);
            [x_srt1(:,t), C_x1] = UKF_FilterStep(x_srt1(:,t), C_x1, newMeasurement, [scale1.mean;[0 0]'],blkdiag(scale1.variance,C_v),l);
            if j == N
                D1 = blkdiag(1,1,0,0,0,0);
                B1 = eye(6) - C_x1*(2*D1*x_srt1(:,t))/((2*D1*x_srt1(:,t))'*C_x1*(2*D1*x_srt1(:,t)))*(2*D1*x_srt1(:,t))';
                x_srt1(:,t) = x_srt1(:,t) + C_x1 * (2*D1*x_srt1(:,t))/((2*D1*x_srt1(:,t))'*C_x1*(2*D1*x_srt1(:,t)))*(1+(2*D1*x_srt1(:,t))'*x_srt1(:,t)-(x_srt1(:,t))'*D1*x_srt1(:,t)-(2*D1*x_srt1(:,t))'*x_srt1(:,t));
                C_x1 = B1 * C_x1 * B1';
                PP(:,:,t) = C_x1;
            end
        end
   %% state prediction
   
        x_srt1(:,t+1) = F_k * x_srt1(:,t);
        C_x1 = F_k * C_x1 * F_k' + C_w;
     end   
     
    x_srt1 = x_srt1(:,1:time_steps);
    [x_fil1new,PP,DD] = urts_smooth1(x_srt1,PP,F_k,C_w);
    
    % calculate the auxiliary function
     [Q(iter),Phh] = c_Q(x_fil1new,PP,Y,C_v,x0,P0,C_w,F_k,DD);
    
     [l] = ccl_1(x_fil1new,Y,C_v,Phh,l);
    if abs(Q(iter) - Q(iter-1)) < tol*abs(Q(iter-1)); break; end  

end

%% illustrate tracking results

for k=1:time_steps
    if mod(k,7) == 0 
        xestnew1 = x_fil1new(:,k);
        h_shape1 = plot_extent([xestnew1(3);xestnew1(5); atan(xestnew1(2)/xestnew1(1));l(1);l(2)], '-','r',.5);
    end
end

%% REM-EOT

J3 = zeros(2);
Nsum = 0;
l = [5 1];
beta = 0.1;
%% initialization
x_srt1 = zeros(n_orientation + 4, time_steps);
x_srt1(:,1) = x0;
C_x1 = P0;

for t = 1:time_steps 
    y = Y{t};
    N = size(y,2);
    if l(2)>l(1)
        l_m = l(2);l(2) = l(1);l(1) = l_m;
    end
    for j = 1 : N
        newMeasurement = y(:,j);
        [x_srt1(:,t), C_x1,l] = UKF_FilterStep(x_srt1(:,t), C_x1, newMeasurement, [scale1.mean;[0 0]'],blkdiag(scale1.variance,C_v),l);
        if j == N
            D1 = blkdiag(1,1,0,0,0,0);
            B1 = eye(6) - C_x1*(2*D1*x_srt1(:,t))/((2*D1*x_srt1(:,t))'*C_x1*(2*D1*x_srt1(:,t)))*(2*D1*x_srt1(:,t))';
            x_srt1(:,t) = x_srt1(:,t)+C_x1*(2*D1*x_srt1(:,t))/((2*D1*x_srt1(:,t))'*C_x1*(2*D1*x_srt1(:,t)))*(1+(2*D1*x_srt1(:,t))'*x_srt1(:,t)-(x_srt1(:,t))'*D1*x_srt1(:,t)-(2*D1*x_srt1(:,t))'*x_srt1(:,t));
            C_x1 = B1 * C_x1 * B1';
            PP(:,:,t) = C_x1;
        end
    end
    
    Phh = HH * C_x1 * HH';
    [l,J3,Nsum] = ccl_2(x_srt1(:,t),Y{t},C_v,Phh,J3,Nsum,l,beta); 
    %% illustrate tracking results

    if mod(t,7) == 0     
        xestnew1 = x_fil1new(:,t);
        h_shape2 = plot_extent([xestnew1(3);xestnew1(5); atan(xestnew1(2)/xestnew1(1));l(1);l(2)], '-.','g',.5);
    end
    %% state prediction
      
    x_srt1(:,t+1) = F_k * x_srt1(:,t);
    C_x1 = F_k * C_x1 * F_k' + C_w;
end

hold on
legend([h_object,  h_shape1, h_shape2], 'Ground Truth', 'EOT-EM', 'EOT-REM')
xlabel('x direction (m)')
ylabel('y direction (m)')
box on
