%Performance for an anisotropic GPS jammer
%Statistics on a large number of Monte-Carlo simulations (50)

clear all

%Definition of R_KF
R_KF=0.08^2; 

%Study of Q_KF
Q_KF=diag([((1))^2 ((1))^2]); 

for l=1:1:1

%Random on P_t (position of the Jammer) - To Uncomment
%Ramdom on x_vec and psi_0 - To Uncomment
%Random on x_t (Power emitted by the Jammer) - To Uncomment

%Initialisation of EKF and UKF
simulation_parameters_anisotropic

P_t_stat(l)=P_t_jammer_num;
x_vec_stat_x(l)=abs(x_vec(1));
x_vec_stat_y(l)=abs(x_vec(2));
psi_0_stat(l)=psi_0;
x_t_vec_stat_x(l)=abs(x_t_vec(1));
x_t_vec_stat_y(l)=abs(x_t_vec(1));

% %Running the file Main-isotropic with EKF
% %Definition of P_cov_ini for EKF
% P_cov_ini=diag([3000^(2) 3000^(2)]);
% Main_anisotropic_extended
% 
% %Recording the important criterion parameters
% if lost_of_simulation == 0
%     total_time_comput_EKF(l)=sum(time_comput_EKF);
%     EKF_fdp(l)=error_RMS_for_CR_extended(N_loops_fb);
%     EKF_final(l)=error_RMS_for_CR_extended(N_loops_vf);
% else
%     total_time_comput_EKF(l)=0;
%     EKF_fdp(l)=0;
%     EKF_final(l)=0;
%     
% end
% %Reinitialising time step and lost of simulation
% k=0; lost_of_simulation = 0;
% 
% %Running the file Main-isotropic with UKF
% %Definition of P_cov_ini for EKF
% P_cov_ini=diag([1200^(2) 1200^(2)]);
% Main_anisotropic_unscented

% Recording the important criterion parameters
% if lost_of_simulation == 0
%     total_time_comput_UKF(l)=sum(time_comput_UKF);
%     UKF_fdp(l)=error_RMS_for_CR_unscented(N_loops_fb);
%     UKF_final(l)=error_RMS_for_CR_unscented(N_loops_vf);
% else
%     total_time_comput_UKF(l)=0;
%     UKF_fdp(l)=0;
%     UKF_final(l)=0;
% end
% 
% Reinitialising time step and lost of simulation
% k=0; lost_of_simulation = 0;
% P_cov_ini=diag([2000^(2) 2000^(2)]);
% 
% defining the number of particles
% N=2000; 
% 
% Running the file particle filter
% Main_anisotropic_particle
% Recording the important criterion parameters
% if lost_of_simulation == 0
%     total_time_comput_particle(l)=sum(time_comput_particle);
%     particle_fdp(l)=error_RMS_for_CR_particle(N_loops_fb);
%     particle_final(l)=error_RMS_for_CR_particle(N_loops_vf);
% else
%     total_time_comput_particle(l)=0;
%     particle_fdp(l)=0;
%     particle_final(l)=0;
% end

%Definition of P_cov_ini for EKF
P_cov_ini=diag([3000^(2) 3000^(2)]);

%Reinitialising time step and lost of simulation
k=0; lost_of_simulation = 0;

%defining the gain for filtering the data
gain=1.1; 

%Running the file particle filter
Main_anisotropic_extended_filtering_data
%Recording the important criterion parameters
if lost_of_simulation == 0
    total_time_comput_EKFF(l)=sum(time_comput_EKFF);
    EKFF_fdp(l)=error_RMS_for_CR_EKFF(N_loops_fb);
    EKFF_final(l)=error_RMS_for_CR_EKFF(N_loops_vf);
else
    total_time_comput_EKFF(l)=0;
    EKFF_fdp(l)=0;
    EKFF_final(l)=0;
end

%Reinitialising time step and lost of simulation
k=0; lost_of_simulation = 0;

%Running the EKF filter with more measurement ratios (technique 1: 10 EKF
%are running at the same time
Main_anisotropic_extended_more_measurement
%Recording the important criterion parameters
if lost_of_simulation == 0
    total_time_comput_EKFPR(l)=sum(time_comput_EKFPR);
    EKFPR_fdp(l)=error_RMS_for_CR_EKF_power_ratio(N_loops_fb);
    EKFPR_final(l)=error_RMS_for_CR_EKF_power_ratio(N_loops_vf);
else
    total_time_comput_EKFPR(l,:)=0;
    EKFPR_fdp(l)=0;
    EKFPR_final(l)=0;
end


%Reinitialising time step and lost of simulation
k=0; lost_of_simulation = 0;
P_r_filt_ratio=0;
%Running the EKF filter with more measurement ratios (technique 1: 10 EKF
%are running at the same time
Main_anisotropic_extended_more_measurement_bis
%Recording the important criterion parameters
if lost_of_simulation == 0
    total_time_comput_EKFPR_bis(l)=sum(time_comput_EKFPR)+sum(time_to_fuse);
    EKFPR_bis_fdp(l)=error_RMS_for_CR_EKF_power_ratio_bis(N_loops_fb);
    EKFPR_bis_final(l)=error_RMS_for_CR_EKF_power_ratio_bis(N_loops_vf);
else
    total_time_comput_EKFPR_bis(l,:)=0;
    EKFPR_bis_fdp(l)=0;
    EKFPR_bis_final(l)=0;
end
end


%%

x_t_vec_mean_x=mean(x_t_vec_stat_x);
x_t_vec_var_x=sqrt(sum(x_t_vec_stat_x.^2)/size(x_t_vec_stat_x,2));


psi_0_mean=mean(psi_0_stat);
psi_0_var=sqrt(sum(psi_0_stat.^2)/size(psi_0_stat,2));


P_t_mean=mean(P_t_stat);
P_t_var=sqrt(sum(P_t_stat.^2)/size(P_t_stat,2));


x_vec_mean_x=mean(x_vec_stat_x);
x_vec_var_x=sqrt(sum(x_vec_stat_x.^2)/size(x_vec_stat_x,2));

% EKF_final_mean=mean(EKF_final);
% EKF_final_var=sqrt(sum(EKF_final.^2)/size(EKF_final,2));
% 
% UKF_final_mean=mean(UKF_final);
% UKF_final_var=sqrt(sum(UKF_final.^2)/size(UKF_final,2));

% particle_final_mean=mean(particle_final);
% particle_final_var=sqrt(sum(particle_final.^2)/size(particle_final,2));

% EKFF_final_mean=mean(EKFF_final);
% EKFF_final_var=sqrt(sum(EKFF_final.^2)/size(EKFF_final,2));
% 
% EKFPR_final_mean=mean(EKFPR_final);
% EKFPR_final_var=sqrt(sum(EKFPR_final.^2)/size(EKFPR_final,2));
% 
% EKFPR_bis_final_mean=mean(EKFPR_bis_final);
% EKFPR_bis_final_var=sqrt(sum(EKFPR_bis_final.^2)/size(EKFPR_bis_final,2));

figure
subplot(2,2,1)
hist(x_t_vec_stat_x)
title('Distribution of the true Position of the GPS Jammer')

subplot(2,2,2)
hist(psi_0_stat)
title('Distribution of the Heading of the GPS Jammer')

subplot(2,2,3)
hist(P_t_stat)
title('Distribution of the Power emitted by the GPS Jammer')

subplot(2,2,4)
hist(x_vec_stat_x)
title('Distribution of the Position of the UAV')

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
% 
% subplot(3,1,1)
% hist(total_time_comput_EKFF)
% title('Distribution of the total simulation time with EKF with filtering on the measure')
% 
% subplot(3,1,2)
% hist(total_time_comput_EKFPR)
% title('Distribution of the total simulation time with EKF with Several Power Ratio - Technique 1')
% 
% subplot(3,1,3)
% hist(total_time_comput_EKFPR_bis)
% title('Distribution of the total simulation time with EKF with Several Power Ratio - Technique 2')

figure
% subplot(3,1,1)
% hist(EKF_final)
% title('Distribution of the Error of the EKF')
% 
% subplot(3,1,2)
% hist(UKF_final)
% title('Distribution of the Error of the UKF')
% 
% subplot(3,1,3)
% hist(particle_final)
% title('Distribution of the Error of the Particle Filter')

subplot(3,1,1)
hist(EKFF_final)
title('Distribution of the Error of the EKF with filtering on the measure')

subplot(3,1,2)
hist(EKFPR_final)
title('Distribution of the Error of the EKF with Several Power Ratio - Technique 1')

subplot(3,1,3)
hist(EKFPR_bis_final)
title('Distribution of the Error of the EKF with Several Power Ratio - Technique 2')