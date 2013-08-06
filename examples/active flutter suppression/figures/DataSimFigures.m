load DataSim
xmax=5;
FS=10;
LW=2;
%% NR=4
figure(1)
subplot(311);
plot(data4.states{1}(:,1), data4.states{1}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with ELLIPSOIDAL approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(312);
plot(data4.states{2}(:,1), data4.states{2}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(313);
plot(data4.states{3}(:,1), data4.states{3}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (PWQ Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');

%%%%%%%%%%%%%%%%%%% Control %%%%%%%%%%%%%%%%%

figure(2)
subplot(311);
plot(data4.ctrl{1}(:,1), data4.ctrl{1}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with ELLIPSOIDAL approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(312);
plot(data4.ctrl{2}(:,1), data4.ctrl{2}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(313);
plot(data4.ctrl{3}(:,1), data4.ctrl{3}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (PWQ Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');



%% NR=6

figure(3)
subplot(311);
plot(data6.states{1}(:,1), data6.states{1}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with ELLIPSOIDAL approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(312);
plot(data6.states{2}(:,1), data6.states{2}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(313);
plot(data6.states{3}(:,1), data6.states{3}(:,2:5),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (PWQ Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');

%%%%%%%%%%%%%%%%%%% Control %%%%%%%%%%%%%%%%%

figure(4)
subplot(311);
plot(data6.ctrl{1}(:,1), data6.ctrl{1}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with ELLIPSOIDAL approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(312);
plot(data6.ctrl{2}(:,1), data6.ctrl{2}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (GLOBAL Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
subplot(313);
plot(data6.ctrl{3}(:,1), data6.ctrl{3}(:,2:3),'LineWidth',LW);
grid on
xlim([0 xmax])
title('BMI approach with QUADRATIC CURVE approximation (PWQ Lyapunov)',...
    'FontSize',FS,...
    'FontName','Arial');
