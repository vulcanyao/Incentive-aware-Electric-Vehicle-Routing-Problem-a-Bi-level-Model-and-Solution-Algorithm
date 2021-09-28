% Result_Visualization

clear


 %%
% cost='200';
% map='c1_#1';
% map='Belgium';
K=3;% Result_Visualization

clear
close all

%%
% startRow = 1
endRow = 300
data_SL = readtable('SL/Sol.csv','Format','%s%f%f%f%f%f%f');
data_1 = readtable('BVRP_gamma_5.0_deltabar_1.5_maxEPI_1/Sol.csv','Format','%s%f%f%f%f%f%f%f%f');
data_2 = readtable('BVRP_gamma_5.0_deltabar_1.0_maxEPI_1/Sol.csv','Format','%s%f%f%f%f%f%f%f%f');
data_3 = readtable('BVRP_gamma_5.0_deltabar_0.5_maxEPI_1/Sol.csv','Format','%s%f%f%f%f%f%f%f%f');
data_4 = readtable('BVRP_gamma_2.5_deltabar_1.5_maxEPI_1/Sol.csv','Format','%s%f%f%f%f%f%f%f%f');
data_5 = readtable('BVRP_gamma_1.5_deltabar_1.5_maxEPI_1/Sol.csv','Format','%s%f%f%f%f%f%f%f%f');



data_SL_obj = table2array(data_SL(1:endRow,3));
data_SL_chargingcost = table2array(data_SL(1:endRow,6));
data_SL_travelcost = table2array(data_SL(1:endRow,7));

data_1_obj = table2array(data_1(1:endRow,3));
data_1_average_deliverysaving = table2array(data_1(1:endRow,7));
data_1_chargingcost = table2array(data_1(1:endRow,8));
data_1_travelcost = table2array(data_1(1:endRow,9));

data_2_obj = table2array(data_2(1:endRow,3));
data_2_average_deliverysaving = table2array(data_2(1:endRow,7));
data_2_chargingcost = table2array(data_2(1:endRow,8));
data_2_travelcost = table2array(data_2(1:endRow,9));

 
data_3_obj = table2array(data_3(1:endRow,3));
data_3_average_deliverysaving = table2array(data_3(1:endRow,7));
data_3_chargingcost = table2array(data_3(1:endRow,8));
data_3_travelcost = table2array(data_3(1:endRow,9));


data_4_obj = table2array(data_4(1:endRow,3));
data_4_average_deliverysaving = table2array(data_4(1:endRow,7));
data_4_chargingcost = table2array(data_4(1:endRow,8));
data_4_travelcost = table2array(data_4(1:endRow,9));


data_5_obj = table2array(data_5(1:endRow,3));
data_5_average_deliverysaving = table2array(data_5(1:endRow,7));
data_5_chargingcost = table2array(data_5(1:endRow,8));
data_5_travelcost = table2array(data_5(1:endRow,9));



% 
 %%  Plot
fontsize = 12
size = 8


figure(1)
for i = 1:(endRow/50)
    data_SL_obj_ele(i) = mean(data_SL_obj(1+50*(i-1) :50*(i)));
    data_1_obj_ele(i) = mean(data_1_obj(1+50*(i-1) :50*(i)));
    data_2_obj_ele(i) = mean(data_2_obj(1+50*(i-1) :50*(i)));
    data_3_obj_ele(i) = mean(data_3_obj(1+50*(i-1) :50*(i)));

    
end
aa = [data_SL_obj_ele',data_1_obj_ele',data_2_obj_ele',data_3_obj_ele'  ];
plot([13:2:23],aa(:,1),'-d','MarkerSize',size,'Color','r','LineWidth',2);
hold on
plot([13:2:23],aa(:,2),'-o','MarkerSize',size,'Color','blue','LineWidth',2);
hold on 
plot([13:2:23],aa(:,3),'-*','MarkerSize',size,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot([13:2:23],aa(:,4),'-*','MarkerSize',size,'Color',[0.6350 0.0780 0.1840],'LineWidth',2);
% set(gca, 'YScale', 'log')
axis([12 24 530 575])
xticks([13:2:23])
legend('Without incentive','With incentive \delta = 1.5 ','With incentive \delta = 1.0 ','With incentive \delta = 0.5 ','Location','northeast')
grid on
xlabel('Size of Network','FontSize',fontsize),ylabel('Average Operation Cost ','FontSize',fontsize)
print('Operation_cost_reduction','-depsc')



figure(2)
for i = 1:(endRow/50)
    data_1_average_deliverysaving_ele(i) = sum(data_1_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_1_average_deliverysaving(1+50*(i-1) :50*(i)));
    data_2_average_deliverysaving_ele(i) = sum(data_2_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_2_average_deliverysaving(1+50*(i-1) :50*(i)));
 data_3_average_deliverysaving_ele(i) = sum(data_3_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_3_average_deliverysaving(1+50*(i-1) :50*(i)));


    
end
aa = [data_1_average_deliverysaving_ele',data_2_average_deliverysaving_ele',data_3_average_deliverysaving_ele'  ]*100/9.05;
bar1 = bar([13:2:23],aa)
bar1(1,3).FaceColor = [0.6350 0.0780 0.1840];
bar1(1,2).FaceColor =[0.4660 0.6740 0.1880];
bar1(1,1).FaceColor = 'blue';
% set(gca, 'YScale', 'log')
% axis([12 20 540 575])
xticks([13:2:23])
legend('With incentive \delta = 1.5 ','With incentive \delta = 1.0 ','With incentive \delta = 0.5 ','Location','northeast')
grid on
xlabel('Size of Network','FontSize',fontsize),ylabel('Average Delivery Fee Saving %','FontSize',fontsize)
print('Delivery_fee_saving','-depsc')







figure(3)
for i = 1:(endRow/50)
    data_SL_obj_ele(i) = mean(data_SL_obj(1+50*(i-1) :50*(i)));
    data_1_obj_ele(i) = mean(data_1_obj(1+50*(i-1) :50*(i)));
    data_4_obj_ele(i) = mean(data_4_obj(1+50*(i-1) :50*(i)));
    data_5_obj_ele(i) = mean(data_5_obj(1+50*(i-1) :50*(i)));

    
end
aa = [data_SL_obj_ele',data_1_obj_ele',data_4_obj_ele',data_5_obj_ele'  ];
plot([13:2:23],aa(:,1),'-d','MarkerSize',size,'Color','r','LineWidth',2);
hold on
plot([13:2:23],aa(:,2),'-o','MarkerSize',size,'Color','blue','LineWidth',2);
hold on 
plot([13:2:23],aa(:,3),'-*','MarkerSize',size,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot([13:2:23],aa(:,4),'-*','MarkerSize',size,'Color',[0.6350 0.0780 0.1840],'LineWidth',2);
% set(gca, 'YScale', 'log')
axis([12 24 525 575])
xticks([13:2:23])
legend('Without incentive','With incentive \gamma = 5 ','With incentive \gamma = 2.5 ','With incentive \gamma = 1.5 ','Location','northeast')
grid on
xlabel('Size of Network','FontSize',fontsize),ylabel('Average Operation Cost ','FontSize',fontsize)
print('Operation_cost_reduction_gamma','-depsc')



figure(4)
for i = 1:(endRow/50)
    data_1_average_deliverysaving_ele(i) = sum(data_1_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_1_average_deliverysaving(1+50*(i-1) :50*(i)));
    data_4_average_deliverysaving_ele(i) = sum(data_4_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_4_average_deliverysaving(1+50*(i-1) :50*(i)));
    data_5_average_deliverysaving_ele(i) = sum(data_5_average_deliverysaving(1+50*(i-1) :50*(i)))/nnz(data_5_average_deliverysaving(1+50*(i-1) :50*(i)));


    
end
aa = [data_1_average_deliverysaving_ele',data_4_average_deliverysaving_ele',data_5_average_deliverysaving_ele'  ]*100/9.05;
bar1 = bar([13:2:23],aa)
bar1(1,3).FaceColor = [0.6350 0.0780 0.1840];
bar1(1,2).FaceColor =[0.4660 0.6740 0.1880];
bar1(1,1).FaceColor = 'blue';
% set(gca, 'YScale', 'log')
% axis([12 20 540 575])
xticks([13:2:23])
legend('With incentive \gamma = 5 ','With incentive \gamma = 2.5 ','With incentive \gamma = 1.5 ','Location','northeast')
grid on
xlabel('Size of Network','FontSize',fontsize),ylabel('Average Delivery Fee Saving %','FontSize',fontsize)
print('Delivery_fee_saving_gamma','-depsc')



























% figure(3)
% for i = 1:(endRow/50)
%     data_SL_chargingcost_ele(i) = mean(data_SL_chargingcost(1+50*(i-1) :50*(i)));
%     data_1_chargingcost_ele(i) = mean(data_1_chargingcost(1+50*(i-1) :50*(i)));
%     
% end
% aa = [data_SL_chargingcost_ele',data_1_chargingcost_ele'  ];
% plot([13:2:23],aa(:,1),'-d','MarkerSize',size,'Color','r','LineWidth',2);
% hold on
% plot([13:2:23],aa(:,2),'-o','MarkerSize',size,'Color','k','LineWidth',2);
% % hold on
% % plot([13:2:17],aa(:,3),'-*','MarkerSize',size,'Color','blue','LineWidth',2);
% % set(gca, 'YScale', 'log')
% % axis([12 20 540 575])
% xticks([13:2:23])
% legend('Without incentive','With incentive ','Location','northeast')
% grid on
% xlabel('Size of Network','FontSize',fontsize),ylabel('Average Charging cost Cost ','FontSize',fontsize)
% print('Charging_cost_saving','-depsc')



