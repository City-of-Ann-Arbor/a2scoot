clf
clc
clear all

%% model COST for each mobility option
mi = [1 1.5 2 2.5 3];

c_umb = [0.01 0.01 0.01 0.01 0.01];      %cost to ride UM bus as a student (free)
c_moped = [0.07 0.10 0.14 0.17 0.21];    %cost calculated by Megan
c_ebike = [0.1 0.15 0.2 0.25 0.31];      %cost calculated by Megan
c_mcycle = [0.29 0.44 0.59 0.74 0.88];   %cost calculated using Megan's car model
c_car = [1.47 2.21 2.95 3.68 4.42];      %cost calculated by Megan
c_aab = [2 2 2 2 2];                     %cost to ride city of AA bus per day ($2/day)
c_scoot = [3.42 4.38 5.98 6.62 6.94];    %base fee + cost/min + fees + taxes
c_uber = [8.09 9.04 10.46 11.73 12.37];  %base fee + cost/min + fees + taxes
c_lyft = [8.42 8.75 9.84 10.88 10.97];   %base fee + cost/min + fees + taxes


one=[c_umb(1) c_moped(1) c_ebike(1) c_mcycle(1) c_aab(1) c_car(1) c_scoot(1) c_lyft(1) c_uber(1)];
one_5=[c_umb(2) c_moped(2) c_ebike(2) c_mcycle(2) c_aab(2) c_car(2) c_scoot(2) c_lyft(2) c_uber(2)];
two=[c_umb(3) c_moped(3) c_ebike(3) c_mcycle(3) c_aab(3) c_car(3) c_scoot(3) c_lyft(3) c_uber(3)];
two_5=[c_umb(4) c_moped(4) c_ebike(4) c_mcycle(4) c_aab(4) c_car(4) c_scoot(4) c_lyft(4) c_uber(4)];
three=[c_umb(5) c_moped(5) c_ebike(5) c_mcycle(5) c_aab(5) c_car(5) c_scoot(5) c_lyft(5) c_uber(5)];


figure(1)
bar([one;one_5;two;two_5;three],'hist')
xlabel('Miles (mi)','fontweight','bold')
ylabel('Cost ($)','fontweight','bold')
legend('UM bus','Moped','e-bike','Motorcycle', 'City bus','Car (gas)', 'Spin e-scooter', 'Lyft', 'Uber','Location','northeastoutside')
set(gca,'XTickLabel',{'1','1.5','2','2.5','3'});


%% model EMISSIONS for each mobility option

e_walk = 37; % gCO2/mi from breathing - not used in graph bc ppl are breathing for every option
e_buson = 52.7;  %emissions released per person riding a hybrid bus during on peak hours = 85(1-.38) gCO2/pmt
e_spinscoot = 0;  %emissions for riding scooter = 0 and Spin offsets emissions for charging/repair
e_rideshare = 244;  %emissions released per person riding uber or lyft (2 passengers minimus so divide avg emissions/pmt driving car by 2
e_car = 367.23;  %avg emissions/pmt driving car
e_truck = 620;  %avg emissions/pmt driving truck
e_busoff = 421.6;  %emissions released per person riding a hybrid bus during off peak hours = 680*(1-.38) gCO2/pmt
e_ebike = 28.148;
e_moped = 88.45;
e_mcycle = 127.89;
e_hybrid = 210.83;
e_EV = 122;
e_scoot = 21.71; %emissions for riding scooter = 0, emissions for charging/repair 86.86 gCO2/passenger-mile

Y = [e_buson e_spinscoot e_rideshare e_car e_truck e_busoff e_ebike e_moped e_mcycle e_hybrid e_EV e_scoot];
X = categorical({'Hybrid Bus (on)','Spin e-Scooter','Rideshare','Personal Car', 'Personal Truck', 'Hybrid Bus (off)', 'e-bike','Moped', 'Motorcycle', 'Hybrid Vehicle','Electric Vehicle', 'e-scooter (other)'});
X = reordercats(X,{'Spin e-Scooter','e-scooter (other)','e-bike','Hybrid Bus (on)','Moped','Electric Vehicle', 'Motorcycle','Rideshare','Hybrid Vehicle','Personal Car', 'Personal Truck', 'Hybrid Bus (off)'});
figure(2)
bar(X,Y)
title('Average U.S. Vehicle Emissions','fontsize',12)
subtitle('Data approximated from literature and may not represent all trip conditions','fontsize',7)
xlabel('Mode of Transportation','fontweight','bold')
ylabel('Emissions (gCO_{2}/PMT)','fontweight','bold')
labels = arrayfun(@(value) num2str(value,'%2.f'),Y,'UniformOutput',false);
text(X,Y,labels,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',8) 

box off

%% Hybrid electric bus emissions / person
