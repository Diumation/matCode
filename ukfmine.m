function [x_state, P_cov, K_UKF_gain] = ukfmine(x_vec_all_1,x_vec_all_k,h_0,P_r_filt_ratio,x_state_ini,P_cov_ini,Q_KF,R_KF)

% Rename the measurement value (for simplicity z)
z=P_r_filt_ratio;                           %Noisy measurement (supposed to be the real one)

% Calculation of the weights for the sigma points selection
L=numel(x_state_ini);                   	%numer of states
m=numel(z);                                 %numer of measurements
alpha=1e-3;                                 %default, tunable
ki=1;                                       %default, tunable  
beta=2;                                     %default, tunable
lambda=alpha^2*(L-ki)-L;                    %scaling factor
c=L+lambda;                                 %scaling factor

Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc=Wm;

Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance

c=sqrt(c);
X=sigmas(x_state_ini,P_cov_ini,c);          %sigma points around x_state_ini

% Unscented transform of the process
[x1,X1,P1,X2]=utprocess(X,Wm,Wc,L,Q_KF);

% Measurement for the sigma points
h=hk(x_vec_all_1,x_vec_all_k,X1,h_0); %Unlinear measurement function

% Unscented transformation of the measurement
[z1,~,P2,Z2]=utmeasure(h,X1,Wm,Wc,m,R_KF);  %unscented transformation of measurments

P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance

K_UKF_gain=P12/(P2);                        %UKF Kalman gain

x_state=x1+K_UKF_gain*(z-z1);               %state update

P_cov=P1-K_UKF_gain*P12';                   %covariance update

function X=sigmas(x_state_ini,P_cov_ini,c)
A = c*chol(P_cov_ini)';
Y = x_state_ini(:,ones(1,numel(x_state_ini)));
X = [x_state_ini Y+A Y-A]; 

function h=hk(x_vec_all_1,x_vec_all_k,x_sigma,h_0)
a=x_sigma(1,:);
b=x_sigma(2,:);
c=x_vec_all_1(1);
d=x_vec_all_1(2);
e=x_vec_all_k(1);
f=x_vec_all_k(2);

n=size(x_sigma,2);
for k=1:n
h(k)=[((a(k)-c)^2+(b(k)-d)^2+(h_0)^2)/((a(k)-e)^2+(b(k)-f)^2+(h_0)^2)];
end

function [y,Y,P_cov_ini,Y1]=utprocess(X,Wm,Wc,n,R)
L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);

for k=1:L                   
    Y(:,k)=X(:,k);
    y=y+Wm(k)*Y(:,k);
end

Y1=Y-y(:,ones(1,L));

P_cov_ini=Y1*diag(Wc)*Y1'+R;

function [y,Y,P_cov_ini,Y1]=utmeasure(f,X,Wm,Wc,n,R)
L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);

for k=1:L
    Y(:,k)=f(:,k);
    y=y+Wm(k)*Y(:,k);
end

Y1=Y-y(:,ones(1,L));

P_cov_ini=Y1*diag(Wc)*Y1'+R;