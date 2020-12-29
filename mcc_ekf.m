close all; clear all; clc;

%% call system dynamic
sysinfo;

%% initiallization
tf = 300;
iter = floor(tf/T); % length of signal
num_of_exprement = 100 ; % number of itterations
num_shot_noise = 15;
start_of_shotnoise = 60;
index_rand_shot = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];

%% Choice the type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = shot noise, 2 = Gaussian mixture noise): ');
    if flag_noise == 1 || flag_noise == 2
        break;
    end
end
    SE_MCC_KF=zeros(num_of_exprement,num_vec,iter);
for Numexper = 1:num_of_exprement
    
   %disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
   
    Q = Q_n1;  R = R_n1;
    MeasErrX = sqrt(Q)*randn(num_vec,iter);
    MeasErrZ = sqrt(R)*randn(num_meas,iter);
    
    % initial conditions for system and estimator(x0 , xhat0)
    x = initial_x; % real states
    xhat2 = initial_x; % (CF): C-Filter 
    xhat3 = initial_x; % (MCC_KF) : Weighted maximum correntropy
    
    % elements of initial estimation error covariance (P0hat)
    P_MCC_KF = initial_P; % (MCC_KF)
    
    
    %% define some arrays to store the main signal and the esitimation signals
    xhat_CF = zeros(num_vec,iter); %(C-Filter)
    xhat_MCC_KF = zeros(num_vec,iter); % (MCC_KF)
    x_main = zeros(num_vec,iter);
    
    %% type of noise
    if flag_noise == 1
      Shot_Noise; 
    else
      Mix_Noise; 
    end
    
    %%
    for t = 1 : 1 : iter
        % make the measurement signal
        z = B*x;
        z = z + MeasErrZ(:,t);

%%  ======================= C_Filter ========================
     xhat2 = A * xhat2  ;
     innov = z - B * xhat2 ;
     norm_innov = sqrt(innov'*innov);
     sigma = 1 * norm_innov;
     K = exp(-(norm_innov^2) /(2*sigma^2));
     xhat2 = xhat2 +  K * B' *(innov) ;
      
     xhat_CF(:,t) = xhat2;
        
%% ======================= run MCC_kf ========================    
    xhat3 = A * xhat3;
    P_MCC_KF = A * P_MCC_KF  * A' + Q;
    invers_R = pinv(R);
    cov_MCC_KF  = P_MCC_KF  ;
    innov = z - B * xhat3;
    norm_innov = sqrt((innov)'*invers_R*(innov));
    sigma = 1 * norm_innov;
    K1 = exp(-(norm_innov^2) /(2*sigma^2));
    Gain = pinv(pinv(cov_MCC_KF )+ K1 * B'*invers_R*B )*K1 * B'*invers_R;
    xhat3 = xhat3 + Gain *(innov);

    P_MCC_KF  = (eye(num_vec) - Gain*B)*cov_MCC_KF *(eye(num_vec) - Gain*B)' + Gain*R*Gain';
    xhat_MCC_KF(:,t)=xhat3;    
    
    %% Simulate the system.
    x = A*x + MeasErrX(:,t);
    x_main(:,t+1) = x;
    end
    x_main(:,iter) = [];
    SE_CF(Numexper,:,:)=(x_main - xhat_CF).^2;
    SE_MCC_KF(Numexper,:,:)=(x_main - xhat_MCC_KF).^2;
    
end

MSE_CF = zeros(num_vec,iter);
MSE_MCC_KF= zeros(num_vec,iter);
%
for i = 1 : iter
    MSE_CF(:,i) = sqrt(mean(SE_CF(:,:,i)))';
    MSE_MCC_KF(:,i) = sqrt(mean(SE_MCC_KF(:,:,i)))';
end
%
MMSE_CF = mean(MSE_CF,2);
MMSE_MCC_KF = mean(MSE_MCC_KF,2);

MMMSE_CF = mean(MMSE_CF);
MMMSE_MCC_KF = mean(MMSE_MCC_KF);

%% Normlized RMSE
max_MCC_KF= max(MSE_MCC_KF.');
max_CF= max(MSE_CF.');
NRMSE = zeros(1,2);

j=1;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_MCC_KF(i,1)/ max_MCC_KF(1,i);
end

j=2;
for i=1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_CF(i,1)/ max_CF(1,i);
end
%% Plot data.
close all
SetPlotOptions;
t = 0 : T : tf-T;
c = 3;
if flag_noise == 1 % ploting in Shot noise case
    figure,
      xlabel('time'), ylabel('Value of shot noise')

    hold on
    X = index_rand_shot;
    Y = MeasErrZ(1,index_rand_shot(1:end));
    for i=1:num_shot_noise
        bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
    end
    
    figure
    hold on
    plot(t,MSE_CF(c,:) ,'r',t,MSE_MCC_KF(c,:) ,'g');legend('C-Filter','MCC-KF'),xlabel('time'), ylabel('Root mean square error')

else % ploting in mixture guassian noise case
      
    f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
    f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));
    
    figure
    subplot(2,1,1),fplot(f1,[-6,6]);title('pdf of process noise')
    subplot(2,1,2),fplot(f2,[-6,6]);title('pdf of measurment noise')
    
    figure
    hold on
    plot(t,MSE_CF(c,:) ,'r',t,MSE_MCC_KF(c,:) ,'g'); legend('C-Filter','MCC-KF'),xlabel('time'), ylabel('Root mean square error')
end
%% result

disp('Mean square error : ');
disp('           x1          x2          x3        x4');
disp(['CF     : ',num2str(MMSE_CF.'),'']);
disp(['MCC_KF : ',num2str(MMSE_MCC_KF.'),'']);