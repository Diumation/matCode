% Comparing the performances of Extended Kalman Filter and Unscented Kalman
% Filter

%% RMSE for EKF and UKF

%For having the same parameters in the Main_isotropic_extended and 
%Main_isotropic_unscented simulations 

%Definition of P_cov_ini
P_cov_ini=diag([4000^(2) 4000^(2)]);

%Definition of Q_KF
Q_KF=diag([((0.1))^2 ((0.1))^2]);  

%Definition of R_KF
R_KF=0.08^2; 

simulation_parameters
%Running the file Main-isotropic with EKF
Main_isotropic_extended
cramer_extended=cramer;
%Reinitialising time step
k=0;
%Running the file Main-isotropic with UKF

%Definition of P_cov_ini
P_cov_ini=diag([1200^(2) 1200^(2)]);

Main_isotropic_unscented
cramer_unscented=cramer;

% Plot the RMSE as a function of the time step
time_step=ind:N_loops_vf;

figure
subplot(3,1,1)
plot(time_step,error_RMS_for_CR_extended(ind:N_loops_vf),'r',time_step,cramer_extended(ind:N_loops_vf),'g');
title('Root Mean Squared Error on the GPS Jammer position - Extended Kalman Filter')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('Extended Kalman Filter','Cramer-Rao Lower Bound')

subplot(3,1,2)
plot(time_step,error_RMS_for_CR_unscented(ind:N_loops_vf),'r',time_step,cramer_unscented(ind:N_loops_vf),'g');
title('Root Mean Squared Error on the GPS Jammer position - Unscented Kalman Filter')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('Unscented Kalman Filter','Cramer-Rao Lower Bound')
subplot(3,1,3)
plot(time_step,error_RMS_for_CR_extended(ind:N_loops_vf),'color',rgb('red')); hold on
plot(time_step,error_RMS_for_CR_unscented(ind:N_loops_vf),'color',rgb('blue')); hold off
title('Root Mean Squared Error on the GPS Jammer position - EKF and UKF')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('Extended Kalman Filter','Unscented Kalman Filter')

%% Influence of P_cov_ini on the EKF and UKF RMSE

clear all

%Definition of Q_KF
Q_KF=diag([((1))^2 ((1))^2]);  

%Definition of R_KF
R_KF=0.08^2; 

%Initialisation of EKF and UKF
simulation_parameters

% Study of P_cov_ini
P_cov_ini_vec(:,:,1)=diag([200^(2) 200^(2)]);
P_cov_ini_vec(:,:,2)=diag([400^(2) 400^(2)]);
P_cov_ini_vec(:,:,3)=diag([600^(2) 600^(2)]); 
P_cov_ini_vec(:,:,4)=diag([800^(2) 800^(2)]); 
P_cov_ini_vec(:,:,5)=diag([1000^(2) 1000^(2)]); 
P_cov_ini_vec(:,:,6)=diag([2000^(2) 2000^(2)]); 
P_cov_ini_vec(:,:,7)=diag([3000^(2) 3000^(2)]); 
P_cov_ini_vec(:,:,8)=diag([4000^(2) 4000^(2)]); 

for l=1:1:8
P_cov_ini=P_cov_ini_vec(:,:,l);

%Running the file Main-isotropic with EKF
Main_isotropic_extended

%Recording the parameters
time_step_extended(l,:)=time_step;
distance_trav_extended(l,:)=distance_trav;
error_RMS_for_CR_extended_vec(l,:)=error_RMS_for_CR_extended;
cramer_extended(l,:)=cramer;

%Reinitialising time step
k=0;
%Running the file Main-isotropic with UKF
Main_isotropic_unscented

%Recording the parameters
time_step_unscented(l,:)=time_step;
distance_trav_unscented(l,:)=distance_trav;
error_RMS_for_CR_unscented_vec(l,:)=error_RMS_for_CR_unscented;
cramer_unscented(l,:)=cramer;

end

% Plotting the RMSE as a function of the time step for the different
% P_cov_ini

time_step=1:N_loops_vf;
figure
plot(time_step_extended(1,time_step),error_RMS_for_CR_extended_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(2,time_step),error_RMS_for_CR_extended_vec(2,time_step),'color',rgb('Blue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(3,time_step),error_RMS_for_CR_extended_vec(3,time_step),'color',rgb('Green'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(4,time_step),error_RMS_for_CR_extended_vec(4,time_step),'color',rgb('Lime'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(5,time_step),error_RMS_for_CR_extended_vec(5,time_step),'color',rgb('Orange'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(6,time_step),error_RMS_for_CR_extended_vec(6,time_step),'color',rgb('Fuchsia'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(7,time_step),error_RMS_for_CR_extended_vec(7,time_step),'color',rgb('Red'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(8,time_step),error_RMS_for_CR_extended_vec(8,time_step),'color',rgb('DarkRed'),'LineStyle','-','LineWidth',2); hold on

plot(time_step_unscented(1,time_step),error_RMS_for_CR_unscented_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(2,time_step),error_RMS_for_CR_unscented_vec(2,time_step),'color',rgb('Blue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(3,time_step),error_RMS_for_CR_unscented_vec(3,time_step),'color',rgb('Green'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(4,time_step),error_RMS_for_CR_unscented_vec(4,time_step),'color',rgb('Lime'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(5,time_step),error_RMS_for_CR_unscented_vec(5,time_step),'color',rgb('Orange'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(6,time_step),error_RMS_for_CR_unscented_vec(6,time_step),'color',rgb('Fuchsia'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(7,time_step),error_RMS_for_CR_unscented_vec(7,time_step),'color',rgb('Red'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(8,time_step),error_RMS_for_CR_unscented_vec(8,time_step),'color',rgb('DarkRed'),'LineStyle',':','LineWidth',2); hold off

title('Root Mean Squared Error on the GPS Jammer position - EKF and UKF')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('\sigma_{ini}= 200','\sigma_{ini}= 400','\sigma_{ini}= 600','\sigma_{ini}= 800',...
    '\sigma_{ini}= 1000','\sigma_{ini}= 2000','\sigma_{ini}= 3000','\sigma_{ini}= 4000')

%RME at the end of the fly-by phase

EKF_fbp(1) = {error_RMS_for_CR_extended_vec(1,N_loops_fb)};
EKF_fbp(2) = {error_RMS_for_CR_extended_vec(2,N_loops_fb)};
EKF_fbp(3) = {error_RMS_for_CR_extended_vec(3,N_loops_fb)};
EKF_fbp(4) = {error_RMS_for_CR_extended_vec(4,N_loops_fb)};
EKF_fbp(5) = {error_RMS_for_CR_extended_vec(5,N_loops_fb)};
EKF_fbp(6) = {error_RMS_for_CR_extended_vec(6,N_loops_fb)};
EKF_fbp(7) = {error_RMS_for_CR_extended_vec(7,N_loops_fb)};
EKF_fbp(8) = {error_RMS_for_CR_extended_vec(8,N_loops_fb)};

disp('RMSE (fly-by phase) - EKF');
disp(EKF_fbp);

UKF_fbp(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_fb)};
UKF_fbp(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_fb)};
UKF_fbp(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_fb)};
UKF_fbp(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_fb)};
UKF_fbp(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_fb)};
UKF_fbp(6) = {error_RMS_for_CR_unscented_vec(6,N_loops_fb)};
UKF_fbp(7) = {error_RMS_for_CR_unscented_vec(7,N_loops_fb)};
UKF_fbp(8) = {error_RMS_for_CR_unscented_vec(8,N_loops_fb)};

disp('RMSE (fly-by phase) - UKF');
disp(UKF_fbp);

%Final RMSE
EKF_final(1) = {error_RMS_for_CR_extended_vec(1,N_loops_vf)};
EKF_final(2) = {error_RMS_for_CR_extended_vec(2,N_loops_vf)};
EKF_final(3) = {error_RMS_for_CR_extended_vec(3,N_loops_vf)};
EKF_final(4) = {error_RMS_for_CR_extended_vec(4,N_loops_vf)};
EKF_final(5) = {error_RMS_for_CR_extended_vec(5,N_loops_vf)};
EKF_final(6) = {error_RMS_for_CR_extended_vec(6,N_loops_vf)};
EKF_final(7) = {error_RMS_for_CR_extended_vec(7,N_loops_vf)};
EKF_final(8) = {error_RMS_for_CR_extended_vec(8,N_loops_vf)};

disp('RMSE (end of simulation) - EKF');
disp(EKF_final);

UKF_final(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_vf)};
UKF_final(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_vf)};
UKF_final(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_vf)};
UKF_final(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_vf)};
UKF_final(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_vf)};
UKF_final(6) = {error_RMS_for_CR_unscented_vec(6,N_loops_vf)};
UKF_final(7) = {error_RMS_for_CR_unscented_vec(7,N_loops_vf)};
UKF_final(8) = {error_RMS_for_CR_unscented_vec(8,N_loops_vf)};

disp('RMSE (end of simulation) - UKF');
disp(UKF_final);

%% Influence of Q_KF on the EKF and UKF RMSE
clear all

%Definition of P_cov_ini
P_cov_ini=diag([400^(2) 400^(2)]);

%Definition of R_KF
R_KF=0.08^2; 

%Keeping a reference for P_t (position of the Jammer)
P_t=0.6

%Keeping a reference for x_t (Power emitted by the Jammer)
x_t_vec=[3.6724e+03    3.4159e+03];

%Keeping a reference for x_vec and psi_0 (position and heading of the UAV)
x_vec=[1.0792e+03 0.8817e+03];

psi_0=1.2510;

%Initialisation of EKF and UKF
simulation_parameters

%Study of Q_KF
Q_KF_vec(:,:,1)=diag([((0.1))^2 ((0.1))^2]); 
Q_KF_vec(:,:,2)=diag([((1))^2 ((1))^2]);  
Q_KF_vec(:,:,3)=diag([((10))^2 ((10))^2]); 
Q_KF_vec(:,:,4)=diag([((100))^2 ((100))^2]); 
Q_KF_vec(:,:,5)=diag([((1000))^2 ((1000))^2]); 

for l=1:1:5
Q_KF=Q_KF_vec(:,:,l)

%Running the file Main-isotropic with EKF
Main_isotropic_extended

%Recording the parameters
time_step_extended(l,:)=time_step;
distance_trav_extended(l,:)=distance_trav;
error_RMS_for_CR_extended_vec(l,:)=error_RMS_for_CR_extended;
cramer_extended(l,:)=cramer;

%Reinitialising time step
k=0;
%Running the file Main-isotropic with UKF
Main_isotropic_unscented

%Recording the parameters
time_step_unscented(l,:)=time_step;
distance_trav_unscented(l,:)=distance_trav;
error_RMS_for_CR_unscented_vec(l,:)=error_RMS_for_CR_unscented;
cramer_unscented(l,:)=cramer;

end


time_step=ind:N_loops_vf;
figure
plot(time_step_extended(1,time_step),error_RMS_for_CR_extended_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(2,time_step),error_RMS_for_CR_extended_vec(2,time_step),'color',rgb('Blue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(3,time_step),error_RMS_for_CR_extended_vec(3,time_step),'color',rgb('Green'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(4,time_step),error_RMS_for_CR_extended_vec(4,time_step),'color',rgb('Lime'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(5,time_step),error_RMS_for_CR_extended_vec(5,time_step),'color',rgb('Orange'),'LineStyle','-','LineWidth',2); hold on

plot(time_step_unscented(1,time_step),error_RMS_for_CR_unscented_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(2,time_step),error_RMS_for_CR_unscented_vec(2,time_step),'color',rgb('Blue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(3,time_step),error_RMS_for_CR_unscented_vec(3,time_step),'color',rgb('Green'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(4,time_step),error_RMS_for_CR_unscented_vec(4,time_step),'color',rgb('Lime'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(5,time_step),error_RMS_for_CR_unscented_vec(5,time_step),'color',rgb('Orange'),'LineStyle',':','LineWidth',2); hold off

title('Root Mean Squared Error on the GPS Jammer position - EKF and UKF')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('\sigma_{state}= 0.1','\sigma_{state}= 1','\sigma_{state}= 10','\sigma_{state}= 100','\sigma_{state}= 1000')

%RME at the end of the fly-by phase

EKF_fbp(1) = {error_RMS_for_CR_extended_vec(1,N_loops_fb)};
EKF_fbp(2) = {error_RMS_for_CR_extended_vec(2,N_loops_fb)};
EKF_fbp(3) = {error_RMS_for_CR_extended_vec(3,N_loops_fb)};
EKF_fbp(4) = {error_RMS_for_CR_extended_vec(4,N_loops_fb)};
EKF_fbp(5) = {error_RMS_for_CR_extended_vec(5,N_loops_fb)};

disp('RMSE (fly-by phase) - EKF');
disp(EKF_fbp);

UKF_fbp(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_fb)};
UKF_fbp(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_fb)};
UKF_fbp(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_fb)};
UKF_fbp(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_fb)};
UKF_fbp(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_fb)};

disp('RMSE (fly-by phase) - UKF');
disp(UKF_fbp);

%Final RMSE
EKF_final(1) = {error_RMS_for_CR_extended_vec(1,N_loops_vf)};
EKF_final(2) = {error_RMS_for_CR_extended_vec(2,N_loops_vf)};
EKF_final(3) = {error_RMS_for_CR_extended_vec(3,N_loops_vf)};
EKF_final(4) = {error_RMS_for_CR_extended_vec(4,N_loops_vf)};
EKF_final(5) = {error_RMS_for_CR_extended_vec(5,N_loops_vf)};

disp('RMSE (end of simulation) - EKF');
disp(EKF_final);

UKF_final(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_vf)};
UKF_final(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_vf)};
UKF_final(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_vf)};
UKF_final(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_vf)};
UKF_final(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_vf)};

disp('RMSE (end of simulation) - UKF');
disp(UKF_final);

%% Influence of R_KF on the EKF and UKF RMSE

clear all

%Definition of P_cov_ini
P_cov_ini=diag([400^(2) 400^(2)]);

%Definition of Q_KF
Q_KF=diag([((1))^2 ((1))^2]);  

%Keeping a reference for P_t (position of the Jammer)
P_t=0.6

%Keeping a reference for x_t (Power emitted by the Jammer)
x_t_vec=[3.6724e+03    3.4159e+03];

%Keeping a reference for x_vec and psi_0 (position and heading of the UAV)
x_vec=[1.0792e+03 0.8817e+03];

psi_0=1.2510;

%Initialisation of EKF and UKF
simulation_parameters

%Study of R_KF
R_KF_vec(:,:,1)=0.01^2; 
R_KF_vec(:,:,2)=0.1^2; 
R_KF_vec(:,:,3)=1^2; 
R_KF_vec(:,:,4)=10^2; 
R_KF_vec(:,:,5)=50^2; 
R_KF_vec(:,:,6)=100^2; 
R_KF_vec(:,:,7)=200^2; 
R_KF_vec(:,:,8)=300^2; 

for l=1:1:8
R_KF=R_KF_vec(:,:,l)

%Running the file Main-isotropic with EKF
Main_isotropic_extended

%Recording the parameters
time_step_extended(l,:)=time_step;
distance_trav_extended(l,:)=distance_trav;
error_RMS_for_CR_extended_vec(l,:)=error_RMS_for_CR_extended;
cramer_extended(l,:)=cramer;

%Reinitialising time step
k=0;
%Running the file Main-isotropic with UKF
Main_isotropic_unscented

%Recording the parameters
time_step_unscented(l,:)=time_step;
distance_trav_unscented(l,:)=distance_trav;
error_RMS_for_CR_unscented_vec(l,:)=error_RMS_for_CR_unscented;
cramer_unscented(l,:)=cramer;

end


time_step=ind:N_loops_vf;
figure
plot(time_step_extended(1,time_step),error_RMS_for_CR_extended_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(2,time_step),error_RMS_for_CR_extended_vec(2,time_step),'color',rgb('Blue'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(3,time_step),error_RMS_for_CR_extended_vec(3,time_step),'color',rgb('Green'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(4,time_step),error_RMS_for_CR_extended_vec(4,time_step),'color',rgb('Lime'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(5,time_step),error_RMS_for_CR_extended_vec(5,time_step),'color',rgb('Orange'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(6,time_step),error_RMS_for_CR_extended_vec(6,time_step),'color',rgb('Fuchsia'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(7,time_step),error_RMS_for_CR_extended_vec(7,time_step),'color',rgb('Red'),'LineStyle','-','LineWidth',2); hold on
plot(time_step_extended(8,time_step),error_RMS_for_CR_extended_vec(8,time_step),'color',rgb('DarkRed'),'LineStyle','-','LineWidth',2); hold on

plot(time_step_unscented(1,time_step),error_RMS_for_CR_unscented_vec(1,time_step),'color',rgb('DarkBlue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(2,time_step),error_RMS_for_CR_unscented_vec(2,time_step),'color',rgb('Blue'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(3,time_step),error_RMS_for_CR_unscented_vec(3,time_step),'color',rgb('Green'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(4,time_step),error_RMS_for_CR_unscented_vec(4,time_step),'color',rgb('Lime'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(5,time_step),error_RMS_for_CR_unscented_vec(5,time_step),'color',rgb('Orange'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(6,time_step),error_RMS_for_CR_unscented_vec(6,time_step),'color',rgb('Fuchsia'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(7,time_step),error_RMS_for_CR_unscented_vec(7,time_step),'color',rgb('Red'),'LineStyle',':','LineWidth',2); hold on
plot(time_step_unscented(8,time_step),error_RMS_for_CR_unscented_vec(8,time_step),'color',rgb('DarkRed'),'LineStyle',':','LineWidth',2); hold off

title('Root Mean Squared Error on the GPS Jammer position - EKF and UKF')
xlabel('Time step')
ylabel('Root Mean Squared Error (m)')
legend('\sigma_{measure}= 0.01','\sigma_{measure}= 0.1','\sigma_{measure}= 1','\sigma_{measure}= 10',...
    '\sigma_{measure}= 50', '\sigma_{measure}= 100','\sigma_{measure}=200','\sigma_{measure}= 300')

%RMSE at the end of the fly-by phase

EKF_fbp(1) = {error_RMS_for_CR_extended_vec(1,N_loops_fb)};
EKF_fbp(2) = {error_RMS_for_CR_extended_vec(2,N_loops_fb)};
EKF_fbp(3) = {error_RMS_for_CR_extended_vec(3,N_loops_fb)};
EKF_fbp(4) = {error_RMS_for_CR_extended_vec(4,N_loops_fb)};
EKF_fbp(5) = {error_RMS_for_CR_extended_vec(5,N_loops_fb)};
EKF_fbp(6) = {error_RMS_for_CR_extended_vec(6,N_loops_fb)};
EKF_fbp(7) = {error_RMS_for_CR_extended_vec(7,N_loops_fb)};
EKF_fbp(8) = {error_RMS_for_CR_extended_vec(8,N_loops_fb)};

disp('RMSE (fly-by phase) - EKF');
disp(EKF_fbp);

UKF_fbp(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_fb)};
UKF_fbp(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_fb)};
UKF_fbp(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_fb)};
UKF_fbp(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_fb)};
UKF_fbp(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_fb)};
UKF_fbp(6) = {error_RMS_for_CR_unscented_vec(6,N_loops_fb)};
UKF_fbp(7) = {error_RMS_for_CR_unscented_vec(7,N_loops_fb)};
UKF_fbp(8) = {error_RMS_for_CR_unscented_vec(8,N_loops_fb)};

disp('RMSE (fly-by phase) - UKF');
disp(UKF_fbp);

%Final RMSE
EKF_final(1) = {error_RMS_for_CR_extended_vec(1,N_loops_vf)};
EKF_final(2) = {error_RMS_for_CR_extended_vec(2,N_loops_vf)};
EKF_final(3) = {error_RMS_for_CR_extended_vec(3,N_loops_vf)};
EKF_final(4) = {error_RMS_for_CR_extended_vec(4,N_loops_vf)};
EKF_final(5) = {error_RMS_for_CR_extended_vec(5,N_loops_vf)};
EKF_final(6) = {error_RMS_for_CR_extended_vec(6,N_loops_vf)};
EKF_final(7) = {error_RMS_for_CR_extended_vec(7,N_loops_vf)};
EKF_final(8) = {error_RMS_for_CR_extended_vec(8,N_loops_vf)};

disp('RMSE (end of simulation) - EKF');
disp(EKF_final);

UKF_final(1) = {error_RMS_for_CR_unscented_vec(1,N_loops_vf)};
UKF_final(2) = {error_RMS_for_CR_unscented_vec(2,N_loops_vf)};
UKF_final(3) = {error_RMS_for_CR_unscented_vec(3,N_loops_vf)};
UKF_final(4) = {error_RMS_for_CR_unscented_vec(4,N_loops_vf)};
UKF_final(5) = {error_RMS_for_CR_unscented_vec(5,N_loops_vf)};
UKF_final(6) = {error_RMS_for_CR_unscented_vec(6,N_loops_vf)};
UKF_final(7) = {error_RMS_for_CR_unscented_vec(7,N_loops_vf)};
UKF_final(8) = {error_RMS_for_CR_unscented_vec(8,N_loops_vf)};

disp('RMSE (end of simulation) - UKF');
disp(UKF_final);

%% Statistics on a large number of Monte-Carlo simulations (50)

clear all

%Definition of R_KF
R_KF=0.08^2; 

%Study of Q_KF
Q_KF=diag([((1))^2 ((1))^2]); 

for l=1:1:150

% File simluation_parameters
%Random on x_t_vec (position of the Jammer) - To Uncomment l 41
%Random on P_t (Power emitted by the Jammer) - To Uncomment l 53
%Ramdom on x_vec and psi_0 - To Uncomment l 71


%Initialisation of EKF and UKF
%Definition of P_cov_ini for EKF
P_cov_ini=diag([4000^(2) 4000^(2)]);

simulation_parameters

P_t_stat(l)=P_t_jammer_num;
x_vec_norm_stat(l)=norm(x_vec);
psi_0_stat(l)=psi_0;
x_t_vec_norm_stat(l)=norm(x_t_vec);

%Running the file Main-isotropic with EKF

Main_isotropic_extended

%Recording the important criterion parameters
total_time_comput_EKF(l)=sum(time_comput_EKF);
EKF_fdp(l)=error_RMS_for_CR_extended(N_loops_fb);
EKF_final(l)=error_RMS_for_CR_extended(N_loops_vf);

%Reinitialising time step
k=0;

%Running the file Main-isotropic with UKF

%Definition of P_cov_ini for UKF
P_cov_ini=diag([1200^(2) 1200^(2)]);
Main_isotropic_unscented

%Recording the important criterion parameters
total_time_comput_UKF(l)=sum(time_comput_UKF);
UKF_fdp(l)=error_RMS_for_CR_unscented(N_loops_fb);
UKF_final(l)=error_RMS_for_CR_unscented(N_loops_vf);

%Reinitialising time step
k=0;
P_cov_ini=diag([2000^(2) 2000^(2)]);
%defining the number of particles

N=2000; 

%Running the file particle filter
Main_isotropic_particle
%Recording the important criterion parameters
total_time_comput_particle(l)=sum(time_comput_particle)+time_generate_particle;
particle_fdp(l)=error_RMS_for_CR_particle(N_loops_fb);
particle_final(l)=error_RMS_for_CR_particle(N_loops_vf);

disp(k)
end

%%

x_t_vec_norm_mean=mean(x_t_vec_norm_stat);
x_t_vec_norm_var=sqrt(sum(x_t_vec_norm_stat.^2)/size(x_t_vec_norm_stat,2));

psi_0_mean=mean(psi_0_stat);
psi_0_var=sqrt(sum(psi_0_stat.^2)/size(psi_0_stat,2));

P_t_mean=mean(P_t_stat);
P_t_var=sqrt(sum(P_t_stat.^2)/size(P_t_stat,2));

x_t_vec_norm_mean=mean(x_t_vec_norm_stat);
x_t_vec_var=sqrt(sum(x_t_vec_norm_stat.^2)/size(x_t_vec_norm_stat,2));

EKF_final_mean=mean(EKF_final);
EKF_final_var=sqrt(sum(EKF_final.^2)/size(EKF_final,2));

UKF_final_mean=mean(UKF_final);
UKF_final_var=sqrt(sum(UKF_final.^2)/size(UKF_final,2));

particle_final_mean=mean(particle_final);
particle_final_var=sqrt(sum(particle_final.^2)/size(particle_final,2));

figure
subplot(2,2,1)
hist(x_t_vec_norm_stat)
title('Distribution of the Position of the GPS Jammer (x-position in m)')

subplot(2,2,2)
hist(psi_0_stat*180/pi)
title('Distribution of the Heading of the GPS Jammer (degrees)')

subplot(2,2,3)
hist(P_t_stat)
title('Distribution of the Power emitted by the GPS Jammer (mW)')

subplot(2,2,4)
hist(x_vec_norm_stat)
title('Distribution of the Position of the UAV (x-position in m)')

% figure
% subplot(3,1,1)
% hist(total_time_comput_EKF)
% title('Distribution of the total simulation time with EKF')
% 
% subplot(3,1,2)
% hist(total_time_comput_UKF)
% title('Distribution of the total simulation time with UKF')
% 
% subplot(3,1,3)
% hist(total_time_comput_particle)
% title('Distribution of the total simulation time with the Particle Filter')

figure
subplot(3,1,1)
hist(EKF_final)
title('Distribution of the Error of the EKF')

subplot(3,1,2)
hist(UKF_final)
title('Distribution of the Error of the UKF')

subplot(3,1,3)
hist(particle_final)
title('Distribution of the Error of the Particle Filter')