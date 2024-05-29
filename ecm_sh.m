clear;
clc;
close all;

% data import
% --------------------------------------------------------------
file_hppc = 'data/HPPC_result.txt';

file_dis_05C = 'data/exp_05c.txt';
file_dis_1C = 'data/exp_1c.txt';
file_dis_2C = 'data/exp_2c.txt';

% 방전 데이터 처리
dis_05C = CellDischargeData.process_discharge_only(file_dis_05C);
dis_1C = CellDischargeData.process_discharge_only(file_dis_1C);
dis_2C = CellDischargeData.process_discharge_only(file_dis_2C);


% Electrical model from HPPC cell data
% --------------------------------------------------------------
data = CellHppcData(file_hppc); %CellHppcData class 호출
params = parameters(); %parameters 호출
ecm = CellEcm(data, params); %CellEcm class 호출

soc = ecm.soc(); 
[v_pts, z_pts] = ecm.ocv(soc, true); %v_pts, z_pts값 반환
coeffs = ecm.curve_fit_coeff(@ecm.func_ttc,5); 
rctau = ecm.rctau_ttc(coeffs);

% Hppc data 검증
soc_hppc = ecm.soc();
ocv_hppc = ecm.ocv(soc_hppc, v_pts, z_pts);

figure;
plot(data.time, ocv_hppc, '.');
ylim([2.6 4.5]);
hold on;

vt_hppc = ecm.vt(soc_hppc, ocv_hppc, rctau);

figure;
plot(data.time, data.voltage, '.', 'DisplayName', 'exp');
hold on;
plot(data.time, vt_hppc, 'DisplayName', 'ecm');
ylim([2.6 4.5]);
legend show;
hold off;

% 방전 데이터 기반 모델 생성
% 0.5C ---------------------------------------
ecm.voltage = dis_05C.voltage;
ecm.time = dis_05C.time;
ecm.current = dis_05C.current;
soc_05C = ecm.soc2(dis_05C);
ocv_05C = ecm.ocv(soc_05C, v_pts, z_pts);
vt_05C = ecm.vt2(dis05C, soc_05C, ocv_05C, rctau);

% 1C -------------------------------------------
ecm.voltage = dis_1C.voltage;
ecm.time = dis_1C.time;
ecm.current = dis_1C.current;
soc_1C = ecm.soc2(dis_1C);
ocv_1C = ecm.ocv(soc_1C, v_pts, z_pts);
vt_1C = ecm.vt2(dis_1C, soc_1C, ocv_1C, rctau);

% 2C -------------------------------------------
ecm.voltage = dis_2C.voltage;
ecm.time = dis_2C.time;
ecm.current = dis_2C.current;
soc_2C = ecm.soc2(dis_2C);
ocv_2C = ecm.ocv(soc_2C, v_pts, z_pts);
vt_2C = ecm.vt2(dis_2C, soc_2C, ocv_2C, rctau);

% plot
figure;
plot(dis_05C.time, dis_05C.voltage, '.', 'DisplayName', 'exp-0.5C');
hold on;
plot(dis_1C.time, dis_1C.voltage, '.', 'DisplayName', 'exp-1C');
plot(dis_2C.time, dis_2C.voltage, '.', 'DisplayName', 'exp-2C');
plot(dis_05C.time, vt_05C, 'DisplayName', 'ecm_0.5C');
plot(dis_1C.time, vt_1C, 'DisplayName', 'ecm_1C');
plot(dis_2C.time, vt_2C, 'DisplayName', 'ecm_2C');
xlabel('Time [s]');
ylabel('Voltage [V]');
legend('Location', 'upper right');
ylim([2.6, 4.5]);
hold off;









