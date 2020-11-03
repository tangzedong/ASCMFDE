clear
close all
radtodeg = 180/pi;
degtorad = pi/180;

BefBAtheta = 30*degtorad;
Votheta = 45*degtorad;
AftBAtheta = 60*degtorad;

BefBAr = [cos(BefBAtheta), sin(BefBAtheta);
    -sin(BefBAtheta), cos(BefBAtheta)];
Vor = [cos(Votheta), sin(Votheta);
    -sin(Votheta), cos(Votheta)];
AftBAr = [cos(AftBAtheta), sin(AftBAtheta);
    -sin(AftBAtheta), cos(AftBAtheta)];

IBefBAr = [cos(-BefBAtheta), sin(-BefBAtheta);
    -sin(-BefBAtheta), cos(-BefBAtheta)];
IVor = [cos(-Votheta), sin(-Votheta);
    -sin(-Votheta), cos(-Votheta)];
IAftBAr = [cos(-AftBAtheta), sin(-AftBAtheta);
    -sin(-AftBAtheta), cos(-AftBAtheta)];

BefBAt = [1,1]';
Vot = [1,1]';
AftBAt = [1,1]';

oBefBAt = IBefBAr*BefBAt;
oVot = IVor*Vot;
oAftBAt = IAftBAr*BefBAt;

T = [oBefBAt, oVot, oAftBAt];
DiffT = oVot - oBefBAt;
subplot(2,3,1)
quiver(zeros(1,3),zeros(1,3), T(1, :),T(2, :), 'Linewidth', 2, 'AutoScale', 'off');
hold on
quiver(oBefBAt(1), oBefBAt(2), DiffT(1), DiffT(2), 'r', 'Linewidth', 2, 'AutoScale', 'off');
hold on
DiffT = BAr*IVor*DiffT;
quiver(oAftBAt(1), oAftBAt(2), DiffT(1), DiffT(2), 'g', 'Linewidth', 2, 'AutoScale', 'off');
xlim([-2,2]);
ylim([-2,2]);
grid on
subplot(2,3,2)
quiver(0,0, BefBAt(1),BefBAt(2), 'Linewidth', 2, 'AutoScale', 'off');
xlim([-2,2]);
ylim([-2,2]);
grid on
subplot(2,3,3)
quiver(0,0, Vot(1),Vot(2), 'Linewidth', 2, 'AutoScale', 'off');
xlim([-2,2]);
ylim([-2,2]);
grid on
subplot(2,3,4)
quiver(0,0, AftBAt(1),AftBAt(2), 'Linewidth', 2, 'AutoScale', 'off');
xlim([-2,2]);
ylim([-2,2]);
grid on