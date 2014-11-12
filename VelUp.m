% Velocity Updates Yo!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
rng(1212);

%Simulate
dt = 0.01; %time differential
t = 0:dt:10;
w = 0.5*pi;
p_gt = sin(w*t);
v_gt = w*cos(w*t);
x_true = [p_gt; v_gt];
l = 1.5;

%Generate measurements
r_meas = p_gt + 1e-1*randn(size(p_gt));
v_meas = v_gt;
meas_idx_v = 1:200:length(p_gt);%find(abs(v_gt) < 1e-3);
%meas_idx_v = [];
meas_idx_r = 1:100:length(p_gt);%find(abs(v_gt) < 1e-3);

r_meas = r_meas(meas_idx_r);
v_meas = v_meas(meas_idx_v);

%Set noise values
sigma_p = sqrt(1e-2);
sigma_v = sqrt(1e-2);
sigma_rm = sqrt(2.5e-3);
sigma_vm = sqrt(2.5e-10);

K = size(x_true,2); %total time steps


%Set variables
A = [1 dt; 0 1];
C_v = [0 1];
C_r = [1 0];
num_r_meas = length(meas_idx_r);
num_v_meas = length(meas_idx_v);


%Set up the matrices
%T
Q = [sigma_p^2 0; 0 sigma_v^2];
R_r = sigma_rm^2;
R_v = sigma_vm^2;
T = sparse(2*K + num_r_meas + num_v_meas, 2*K + num_r_meas + num_v_meas);
for qi = 1:K
    T((qi*2 -1):qi*2, (qi*2 -1):qi*2) = Q;
end
for ri = 1:(num_r_meas)
    T((ri +2*K), (ri+ 2*K)) = R_r;
end
for vi = 1:(num_v_meas)
    T((vi + 2*K + num_r_meas), (vi + 2*K + num_r_meas)) = R_v;
end



%G
i_idx = 1:2*K;
j_idx = 1:2*K;
G_values = zeros(1,2*K);
i_idx = [i_idx 2*K+1:2*K+(num_r_meas + num_v_meas)];
j_idx = [j_idx 2*K+1:2*K+(num_r_meas + num_v_meas)];
G_values = [G_values ones(1,num_r_meas + num_v_meas)];
G  = sparse(i_idx,j_idx,G_values);

%H
H = sparse(2*K + num_r_meas + num_v_meas, 2*(K + 1));
for ai = 1:K
    H((ai*2 -1):(ai*2), (ai*2 -1):ai*2+2) = [-A eye(2)];
end

for ri = 1:(num_r_meas)
    H(ri+2*K, 2*meas_idx_r(ri)+1:2*meas_idx_r(ri)+2) = C_r;
end
for vi = 1:(num_v_meas)
    H(vi+2*K + num_r_meas, (2*meas_idx_v(vi)+1):(2*meas_idx_v(vi)+2)) = C_v;
end

% Compute
% Form matrices 
y_obs = [r_meas'; v_meas'];
y_obs = y_obs(:);
z = [zeros(2*K,1); y_obs];

Amat = H'/T*H;
bvec = H'/T*G*z;

%Solve!
x_star =  Amat\bvec; %Fast!
x_star_var = diag(inv(Amat)); %Slow, but only necessary for uncertainty

% Plot
% Extract position and velocity from k=2 to end
p_star = x_star(3:2:end);
v_star = x_star(4:2:end);
p_star_var = x_star_var(3:2:end);
v_star_var = x_star_var(4:2:end);


close all;
figure
p_errors = p_star' - x_true(1,:);
%plot(t, x_errors, 'b');
subplot(2,1,1)
plot(t, p_star', 'b');
hold on;
plot(t, x_true(1,:), '--k');
plot(t, p_star - 3*sqrt(p_star_var), '-r', 'Linewidth', 1);
plot(t, p_star + 3*sqrt(p_star_var), '-r', 'Linewidth', 1);

xlabel('Time [s]');
ylabel('Position [m]');
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12)

subplot(2,1,2)
plot(t, v_star', 'b');
hold on;
plot(t, x_true(2,:), '--k');
plot(t, v_star - 3*sqrt(v_star_var), '-r', 'Linewidth', 1);
plot(t, v_star + 3*sqrt(v_star_var), '-r', 'Linewidth', 1);

xlabel('Time [s]');
ylabel('Speed [m/s]');
set(gca,'FontSize',12);
set(findall(gcf,'type','text'),'FontSize',12)

% Histogram
figure
hist(p_errors)
ylabel('Frequency as Occurences');
xlabel('Position Error [m]');
title(sprintf('Position Error \n 3-sigma: %.5f [m]', 3*std(p_errors))); 
set(gca,'FontSize',12)
set(findall(gcf,'type','text'),'FontSize',12)
pbaspect([2 1 1])
