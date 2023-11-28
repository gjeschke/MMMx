data2 = load('MVN2_flexibility_chain_A.csv');
data4 = load('MVN4_flexibility_chain_A.csv');

figure(1); clf; hold on;
h2 = plot(data2(:,1),data2(:,2),'.','MarkerSize',14);
h4 = plot(data4(:,1),data4(:,2),'.','MarkerSize',14);
set(gca,'FontSize',12);
axis([0,max(data2(:,1))+1,0,1]);
legend([h2,h4],'e002','e004','Location','northwest');
xlabel('Residue');
ylabel('Site-specific flexibility f_i');

order2 = load('MVN2_order_chain_A.csv');
order4 = load('MVN4_order_chain_A.csv');

figure(2); clf; hold on;
h2 = plot(order2(:,1),order2(:,2),'.','MarkerSize',14);
h4 = plot(order4(:,1),order4(:,2),'.','MarkerSize',14);
set(gca,'FontSize',12);
axis([0,max(data2(:,1))+1,0.19,0.225]);
legend([h2,h4],'e002','e004','Location','northwest');
xlabel('Residue');
ylabel('Site-specific order o_i');
