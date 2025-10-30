function [x_state, P_cov, K_EKF_gain] = KF_form(x_vec_all_1,x_vec_all_k,h_0,P_r_filt_ratio,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF)

% Extended Kalman Filter

%% Prediction

%% Equation 1 - Initialisation
x_state_s = F_KF*x_state_ini;

%% Equation 2
P_cov_s = F_KF*P_cov_ini*transpose(F_KF)+G_KF*Q_KF*transpose(G_KF); % Error covariance extrapolation

%% Correction

%% nonlinear measurement eq
h = hk(x_vec_all_1,x_vec_all_k,x_state_s,h_0);

%% Jacobian of nonlinear measurement eq.
H = JH(x_vec_all_1,x_vec_all_k,x_state_s,h_0);

%% Equation 3: Innovation
v = (P_r_filt_ratio-h);

%% Equation 4: Innovation
S = H*P_cov_s*transpose(H)+R_KF;

%% Equation 5: Kalman gain
K_EKF_gain = P_cov_s*transpose(H)*(S)^-1;

%% Equation 6: State update
x_state = x_state_s + K_EKF_gain*v;

%% Equation 7: Error covariance update
P_cov = (eye(2)-K_EKF_gain*H)*P_cov_s;

%% h(X): Nonlinear measurement eq
function h=hk(x_vec_all_1,x_vec_all_k,x_state_s,h_0)
a=x_state_s(1);
b=x_state_s(2);
c=x_vec_all_1(1);
d=x_vec_all_1(2);
e=x_vec_all_k(1);
f=x_vec_all_k(2);

h=[((a-c)^2+(b-d)^2+(h_0)^2)/((a-e)^2+(b-f)^2+(h_0)^2)];


%% partial_h/partial_X: Jacobian of measurement eq
function H=JH(x_vec_all_1,x_vec_all_k,x_state_s,h_0)
a=x_state_s(1);
b=x_state_s(2);
c=x_vec_all_1(1);
d=x_vec_all_1(2);
e=x_vec_all_k(1);
f=x_vec_all_k(2);

num=(a-c)^2+(b-d)^2+(h_0)^2;
den=(a-e)^2+(b-f)^2+h_0^2;

H=[(2*(a-c)*den-2*(a-e)*num)/den^2 (2*(b-d)*den - 2*(b-f)*num)/den^2];


