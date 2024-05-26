# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import params

import cell_discharge_data
import cell_temperature_data
import cell_hppc_data
import cell_ecm
import thermal_model
import utils

# Data files
# ----------------------------------------------------------------------------
file_hppc = 'data/HPPC_result.txt'

file_dis_05c = 'data/exp_05c.txt'
file_dis_1c = 'data/exp_1c.txt'
file_dis_2c = 'data/exp_2c.txt'

# file_temp_1c = 'data/cell-discharge-temperature-1c.lvm'
# file_temp_2c = 'data/cell-discharge-temperature-2c.lvm'
# file_temp_3c = 'data/cell-discharge-temperature-3c.lvm'

# Processed cell discharge data for just the discharge section
# ----------------------------------------------------------------------------

dis_05c = cell_discharge_data.CellDischargeData.process_discharge_only(file_dis_05c)
dis_1c = cell_discharge_data.CellDischargeData.process_discharge_only(file_dis_1c)
dis_2c = cell_discharge_data.CellDischargeData.process_discharge_only(file_dis_2c)





# temp_1c = cell_temperature_data.CellTemperatureData.process(file_temp_1c, dis_1c.ti, dis_1c.tf)
# temp_2c = cell_temperature_data.CellTemperatureData.process(file_temp_2c, dis_2c.ti, dis_2c.tf)
# temp_3c = cell_temperature_data.CellTemperatureData.process(file_temp_3c, dis_3c.ti, dis_3c.tf)

# Electrical model from HPPC cell data
# ----------------------------------------------------------------------------
data = cell_hppc_data.CellHppcData(file_hppc) # hppc 데이터 불러와서 data 에 저장
ecm = cell_ecm.CellEcm(data, params)  # hppc 데이터 기반으로 ecm ()
soc = ecm.soc()
_, _, _, v_pts, z_pts = ecm.ocv(soc, pts=True)
coeffs = ecm.curve_fit_coeff(ecm.func_ttc, 5)
rctau = ecm.rctau_ttc(coeffs)



# HPPC data validation
soc_hppc = ecm.soc()
ocv_hppc = ecm.ocv(soc_hppc, vz_pts=(v_pts, z_pts))
fig, ax = plt.subplots()
ax.plot(data.time, ocv_hppc, marker='.', label='ocv')
plt.ylim([2.6, 4.5])  

vt_hppc = ecm.vt(soc_hppc, ocv_hppc, rctau)

fig, ax = plt.subplots()
ax.plot(data.time, data.voltage, marker='.', label='exp')
ax.plot(data.time, vt_hppc, label='ecm')
plt.ylim([2.6, 4.5])     # Y축의 범위: [ymin, ymax]
plt.show()

# Thermal model from Discharge 0.05C
# ----------------------------------------------------------------------------

# ocv_input=pd.read_csv('data/ocv.csv')
# ecm.time =np.array(ocv_input).T[0]
# ecm.current =np.array(ocv_input).T[1]
# soc_005c = ecm.soc()
# ocv_005c = ecm.ocv(soc_005c, vz_pts=(v_pts, z_pts))
# vt_005c = ecm.vt(soc_005c, ocv_005c, rctau)

# fig, ax = plt.subplots()
# ax.plot(ecm.time, vt_005c, marker='.', label='voltage')
# ax.plot(ecm.time, ocv_005c, label='ocv')
# utils.config_ax(ax, xylabels=('Time [s]', 'Voltage [V]'), loc='upper right')



#Thermal model from Discharge 05C
#----------------------------------------------------------------------------
ecm.voltage = dis_05c.voltage
ecm.time = dis_05c.time
ecm.current =dis_05c.current
soc_05c = ecm.soc2(dis_05c)
ocv_05c = ecm.ocv(soc_05c, vz_pts=(v_pts, z_pts))
vt_05c= ecm.vt2(dis_05c,soc_05c, ocv_05c, rctau)

#Thermal model from Discharge 1C
#----------------------------------------------------------------------------
ecm.voltage = dis_1c.voltage
ecm.time = dis_1c.time
ecm.current =dis_1c.current
soc_1c = ecm.soc2(dis_1c)
ocv_1c = ecm.ocv(soc_1c, vz_pts=(v_pts, z_pts))
vt_1c= ecm.vt2(dis_1c,soc_1c, ocv_1c, rctau)


#Thermal model from Discharge 2C
#----------------------------------------------------------------------------
ecm.voltage = dis_2c.voltage
ecm.time = dis_2c.time
ecm.current =dis_2c.current
soc_2c = ecm.soc2(dis_2c)
ocv_2c = ecm.ocv(soc_2c, vz_pts=(v_pts, z_pts))
vt_2c= ecm.vt2(dis_2c, soc_2c, ocv_2c, rctau)


# Plot
# ----------------------------------------------------------------------------

fig, ax = plt.subplots()
ax.plot(dis_05c.time, dis_05c.voltage, marker='.', label='exp-0.5c')
ax.plot(dis_1c.time, dis_1c.voltage, marker='.', label='exp-1c')
ax.plot(dis_2c.time, dis_2c.voltage, marker='.', label='exp-2c')
ax.plot(dis_05c.time, vt_05c, label='ecm_0.5c')
ax.plot(dis_1c.time, vt_1c, label='ecm_1c')
ax.plot(dis_2c.time, vt_2c, label='ecm_2c')
utils.config_ax(ax, xylabels=('Time [s]', 'Voltage [V]'), loc='upper right')
plt.ylim([2.6, 4.5]) 

plt.show()