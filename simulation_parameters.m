%   -----------------------------------------------------------------------
%   --------------------    Simulation parameters    ----------------------
%   -----------------------------------------------------------------------
%   Parameters (Simulation design): Change according to desired simulation
%   Within the simulation, these are fixed

%   Area parameters and number of iterations are global variables
    global x_bnd y_bnd N_loops_fb
    
    global d2r

%   Constants (not to be modified)
d2r=pi/180;                                                                 %   Value in rad = Value in deg * d2r
r2d=1/d2r;                                                                  %   Value in deg = Value in rad * r2d
g_0=9.81;                                                                   %   Gravity acceleration assumd constant [m/s^2]
f_L1=1575.42*10^6;                                                          %   L1 frequency [Hz]
c_0=299792458;                                                              %   speed of light [m/s]
k_b=1.3806488*(10^(-23));                                                   %   Boltzmann constant

    %   Simulation parameters
        D_T=1/2;                    %   Sampling time of the simulation (for Euler integration) [s] - half a second is fine for filter accuracy and animation speed
        t_0=0;                                                              %   Initial time [s]
        t_f_fb=4.5*60;              %   Final time for flyby [s] - Check with speed to know what distance will be travelled - 6 minutes at 30m/s is fine for straight flyby
        t=(t_0:D_T:t_f_fb-t_0)';                                         	%   Time vector for flyby [s]
        N_loops_fb=size(t,1);                                              	%   Number of loops in the flyby simulation
        t_f_vf=15*60;               %   Final time for the Vector field part. Note it is not its duration but the time at which it stops
        t=(t_0:D_T:t_f_vf-t_0)';                                            %   Time vector for whole simulation
        N_loops_vf=size(t,1);                                               %   Total number of loops for the simulation
        
        
        
    %   Search Area parameters
        x_bnd=12*10^3;                                                      %   x area boundary [m]
        y_bnd=12*10^3;                                                      %   y area boundary [m]
        A_area=x_bnd*y_bnd;                                                 %   Area of search [m²]
    
        
    %   Jammer parameters
   	%   Static Jammer true location : located within a square centred
   	%   inside the search area. These parameters are not known by the UAV
        x_t_vec=place_jammer();                                             %   See corresponding function. It places the jammer randomly in a square in the search area
        %   GPS jammer model for the simulation
         P_t_min=1*(10^(-3));                                                %   [W] - Generally around 1mW
         P_t_max=650*(10^(-3));                                              %   [W] - Generally 650mW
        %   Jammer power in the L1 band: to be adjusted for desired
        %   simulation. Assumed always constant (civil jammers)
        
        %   Jammer orientation
        psi_jammer=0;                                                      %   Degrees [0-360]
        psi_jammer=psi_jammer*d2r;                                          %   In radians
        
        %   Jammer power             
        P_t=(P_t_min+(P_t_max-P_t_min)*rand(1,1));                          %   Power [W] - Line to comment if the user wants to specify a given power within power bounds
        %   Display jammer power as it is an important simulation parameter
        P_t_jammer_num=P_t*(10^3);
        P_t_jammer_str=num2str(P_t*(10^3));                                 %   convert to mW, then to string for display
        disp('--> GPS jammer Power in L1 band for this simulation is :');   %   display
        disp([P_t_jammer_str,' mW']);                                       %   display
        
        %   Jammer gain
        G_t=1;                                                              %   Isotropic, lossless antenna []
        %   Jammer sweep range: [f_L1-f_min ; f_L1+f_max]
        f_min=f_L1-30*10^6;                                                 %   [Hz] - f_min and f_max around 10-20-30 MHz
        f_max=f_L1+30*10^6;                                                 %   [Hz]
    
   
        
        
    %   UAV parameters    
        %   Initial position and heading
        [x_vec psi_0]=place_uav();	%   See corresponding function. It places the UAV randomly in a small square in the South-West area with a random heading                          
        
        
        %   Altitude-hold
    	h_0=125;                                                            %   Constant altitude of the UAV [m]        
        
        %   UAV airspeed
        V_g=28.3;                                                           %   Average Cruising Airspeed of the aerosonde [m/s]
        V_min=15;                                                           %   Minimum safe speed [m/s]
        V_max=50;                                                           %   Maximum safe speed [m/s]
        %   UAV minimum turn radius
        min_turn_r=200;                                                       %   Minimum safe turn radius [m]
        %   UAV Antenna Gain
        G_r=1;                                                              %   Isotropic, lossless antenna []
        %   Obtain side on which the jammer lies: 
        %   port of UAV--> true_side=1     starboard of UAV--> true_side=0
        true_side=get_true_side(x_t_vec,x_vec,psi_0);                       %   See corresponding function for details
        
        
    %   UAV path parameters: the winding path is generated by cosinusoidal
    %   heading command   (psi-psi_0)=psi_range*cos((2pi/dist_period)*distance_travelled)
        dist_period=(10/4)*10^3;                                                 %   Distance travelled by UAV during one winding path period [m] 
        psi_range=(2*pi/4);                                                 %   Psi variation range during winding path [rad]
    
        
    %   Measurement process parameters
        %   Friis' equation constant parameter called gamma_0
        gamma_0=G_t*G_r*((c_0/(4*pi*f_L1))^2);                              %   Coefficient assumed constant
        %   Measurement noise
        Temperature=23;                                                     %   Temperature of the sensor [C°]
        P_thermal_noise=k_b*(Temperature+273.15)*abs(f_max-f_min);          %   Thermal noise using Johnson–Nyquist equation
        P_thermal_noise_dBm=10*log10(1000*P_thermal_noise);                 %   Converstion in dBm
        %   Filtering
        low_pass_freq=0.06;                                                 %   Low pass filter cut-off frequency W_n: check help butter for more information (good values 0.01 - 0.1)
        butter_order=2;                                                     %   Order of the low-pass filter generated by the butterworth command
        [b_butter,a_butter]=butter(butter_order,low_pass_freq);             %   Obtain Butterworth filter coefficients
        
        %   Sensor power measurement standard deviation: Choose and disable
        %   lines for noise specification or for thermal noise approximation:
        sig_P_r=-120;                                                       %   Gaussian noise affects Power measurement P_r: std sigma in [dBm]
        %sig_P_r=3*P_thermal_noise_dBm;                                     %   std sigma based on thermal noise [dBm]
        %sig_P_r_W=P_thermal_noise/3;                                       %   std sigma based on thermal noise [W]
        sig_P_r_W=((10^(sig_P_r/10))/1000);                                 %   std sigma in [W]
        
        
    %   Geolocation process parameters
        %   Confidence interval/max band for received power peak determination
            conf_intvl=(1-(2/100));                                         %   Default setting: 2%
            k_H_g=0;                                                        %   Serves as a check condition on the passing of a minimum distance to the jammer for straight flyby
            P_r_filt_max=0;                                                 %   Initialise maximum filtered power received. Serves as check condition for max determination for straight flyby
            band_found=0;                                                   %   Serves as boolean to check whether the max band has been found straight flyby
            
    	%   Extended Kalman filtering based on alpha-measurements
        %   Parameters for EKF (UKF) start
            k_C_prim=floor((0.5/100)*N_loops_fb)+1;                         %   Initialise k_C_prim: step used to determine whether two circles seperated by a given distance intersect
            obs_check=0;                                                    %   Boolean parameter for EKF (UKF) start: initialise start as false
            k_obs=1;                                                        %   Step at which obs_check turns true
            div_EKF_bool=0;                                                 %   Boolean for EKF (UKF) track divergence: initialise as false
        %   EKF (UKF) initialisation:
            F_KF=eye(2);                                                    %   Dynamics matrix: unity because model is static 
            G_KF=eye(1);                                                    %   Noise matrix: unity for pure additive gaussian noise
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---------- < Q_KF must be set up appropriately > ------------%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Q_KF=diag([((1))^2 ((1))^2]);                             %   Process noise matrix: better to be small std for position and power
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---------- < Q_KF must be set up appropriately > ------------%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R_KF=0.08^2;                                                	%   Specify noise on alpha: enable if wanted  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            x_state_ini=[x_bnd/2 1*y_bnd/2]';                                 %   Initial state guess - Middle of the area is the first guess
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---------- < Q_KF must be set up appropriately > ------------%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P_cov_ini=diag([400^(2) 400^(2)]);                            %   Initial state covariance guess - Change if needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        %    Parameters for EKF (UKF) rerun
            re_run_bool=0;     	
            
    %   Vector field part parameters
        %   Desired radius r_d for the orbit around the jammer
        r_d=1000;                                                           %   in metres [m]
        alpha_f=1;                                                          %   Normalising constant on V_g during orbital adaptation: alpha=1 --> V=V_g
        K_LVFG_psi=3.0;                                                     %   Gain on the heading rate command for the LVFG      
            
 	%   Prepare plot for animation
    global plot_scaling fig_offset 

        %   Plot parameters
        N_plots=1;              %   Counter for plots
        plot_scaling=10^3;      %   IMPORTANT: Scale to kilometers --> Changing this parameter will need changing axes labels
        fig_offset=1/100;       %   Area offset to figure boundaries [% of bounds x_bnd & y_bnd]
        p_e=95/100;             %   Confidence probability for the error covariance based on P_cov
        
%%

%   -----------------------------------------------------------------------
%   ------    Initialise simulation parameters for flyby phase   ----------
%   -----------------------------------------------------------------------


%   UAV
x_vec_all=zeros(N_loops_fb,2);                                                 %   Initialise uav position vector
x_vec_all(1,:)=x_vec;                                                          %   idem for initial postion

x_vec_dot=zeros(N_loops_fb,2);                                                 %   Initialise uav velocity vector: derivative of x_vec_all

psi_all=zeros(N_loops_fb,1);                                                   %   Initialise uav heading vector
psi_all(1,:)=psi_0;                                                            %   Idem for initial heading
psi_dot=zeros(N_loops_fb,1);                                                   %   Initialise uav heading derivative vector

jammer_UAV_vec_p=zeros(N_loops_fb,3);                                          %    Initialise Jammer-->UAV vector (3D)
elev_angle=zeros(N_loops_fb,1);                                                %    Initialise elevation angle to vertical
azimuth_angle=zeros(N_loops_fb,1);                                             %    Initialise jammer-->UAV azimuth angle to horizontal
azimuth_rel_angle=zeros(N_loops_fb,1);                                         %    Initialise jammer-->UAV azimuth angle to jammer direction


%   Measurements
r_true=zeros(N_loops_fb,1);                                                    %   True slant range
P_r_true=zeros(N_loops_fb,1);                                                  %   True received power
P_r_meas=zeros(N_loops_fb,1);                                                  %   Measured received power: with noise

%   Processing
r_est_l=zeros(N_loops_fb,1);                                                   %   Lower range estimation
r_est_h=zeros(N_loops_fb,1);                                                   %   Upper range estimation
r_est_unc=zeros(N_loops_fb,1);                                                 %   Range estimation uncertainty
P_r_filt_ratio=zeros(N_loops_fb,1);                                            %   Alpha: Power ratio between initial and current : see 'alpha' in report
centre_geo_circle=zeros(N_loops_fb,2);                                         %   Centre of geolocation circle at instant k
radius_geo_circle=zeros(N_loops_fb,1);                                         %   Radius of geolocation circle at instant k

%   Filters
    %   EKF or UKF
    x_state=zeros(2,N_loops_fb);                                             	%   Updated EKF (UKF) state vector for all steps                               
    P_cov=zeros(2,2,N_loops_fb);                                             	%   EKF (UKF) Covariance matrix for all
    K_EKF_gain=zeros(2,N_loops_fb);                                           	%   Kalman gain storage

    
                                        
%   Simulation data
d_uav=zeros(N_loops_fb,1);                                                     %   Distance travelled by the UAV

