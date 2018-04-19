medi_store = load('lambdas_36_midn.mat');
medi_store = medi_store.medi_store;
price_store = load('price_36_midn.mat');
price_store = price_store.price_store;
n_period = 24;
problem = cell(1,n_period);
dp_result = cell(1,n_period);
dp_time = zeros(1,n_period);
adp_result = cell(1,n_period);
adp_time = zeros(1,n_period);
adp_firstD = zeros(n_period,1);
time_slot = 12;
arrive_time =18;
departure_time = 9;
post_decision_dp =ones(n_period,3);
pre_decision_dp = ones(n_period+1,1);
post_decision_adp =ones(n_period,1);
pre_decision_adp = ones(n_period+1,1);
pre_decision_adp_rounded = ones(n_period+1,1);
pre_decision_dp(1,:) = [75];
pre_decision_adp(1,:) = [77];
pre_decision_adp_rounded(1,:) = [77];
forecast = medi_store;
%forecast = medi_store + mod(randn(36,3),3);
%forecast(:,1) = Generate_ARMA_Forecast_HEMS(medi_store(:,1),0.05,36);
%forecast(:,2) = Generate_ARMA_Forecast_HEMS(medi_store(:,2),0.1,36);
%forecast(:,3) = round(Generate_ARMA_Forecast_HEMS(medi_store(:,3),0.1,36));
%forecast(:,3) = forecast(:,3)-mod(forecast(:,3),4);
%forecast(:,2) = round(forecast(:,2));
H_arrive_time = arrive_time+1;
H_departure_time = departure_time+1;
day = 2;
hour = 0;

for i = 1:n_period
    medi = medi_store(i:i+time_slot-1,:);
    price = price_store(:,i:i+time_slot-1);
    H_arrive_time = H_arrive_time - 1;
    H_departure_time = H_departure_time - 1;
    if H_departure_time == 0
        H_departure_time = 24;
    end
    ad_time = [H_arrive_time,H_departure_time];
    
    %ADP
    [problem,result,time] = HEMS_Main(false,time_slot,pre_decision_adp_rounded(i,:),medi,price,ad_time,151,60,30);
    adp_firstD(i,:) = result.first_decision;
    
    %DP
    %To run DP, Comment above two lines and uncomment below two lines
%    [problem,result,time] = HEMS_Main(true,time_slot,pre_decision_adp_rounded(i,:),medi,price,ad_time);
%    adp_firstD(i,:) = cell2mat(result.dpbi_policy(1));
    
    s_range_mat = cell2mat(problem.params.appliance_range);
    adp_time(1,i) = time;
    post_decision_adp(i,:) = pre_decision_adp(i,:)+adp_firstD(i,:);
    pre_decision_adp(i+1,1) = 0.9*pre_decision_adp(i,1)+ 0.1*forecast(i,1)+adp_firstD(i,1);
    pre_decision_adp_rounded(i+1,1) = min(round(pre_decision_adp(i+1,1),1),s_range_mat(1,2));
end

combined_adp = ones(n_period*2+1,3);
for n = 1:n_period
    combined_adp(2*n-1,:) = pre_decision_adp(n,:);
    combined_adp(2*n,:) = post_decision_adp(n,:);
    x(2*n-1) = n;
    x(2*n) = n+0.2;
    x_r(n) = n+0.6;
end
combined_adp(2*n_period+1,:) = pre_decision_adp(n_period+1,:);
x(2*n_period+1) = n_period+1;


yyaxis left
rectangle('Position',[14-1 0 6 85],'EdgeColor',[1 1 1],'FaceColor',[1,0,1]);
hold on;
rectangle('Position',[11-1 0 3 85],'EdgeColor',[1 1 1],'FaceColor',[1,0.8,1]);
hold on;
rectangle('Position',[20-1 0 2 85],'EdgeColor',[1 1 1],'FaceColor',[1,0.8,1]);
hold on;
plot(x-1,combined_adp(:,1),'black','DisplayName','Room Temperature','LineWidth',2);
hold on;
%plot(x,combined_dp(:,1),'DisplayName','dp_Room Temperature');
hold on;
line2 = refline([0 75]);
line2.LineStyle = '-.';
line2.Color = 'g';
line2.DisplayName = 'Desired Temperature';
ylabel('Room and Desired Temperature(F)');
ylim([60 85]);
hold on;
yyaxis right
bar(x_r-1,forecast(1:24,1),'DisplayName','Outside Temperature');
ylabel('Outside Temperature(F)');
ylim([75 130]);
xlabel('Time(hr)');
xlim([0 24]);
legend('show');
