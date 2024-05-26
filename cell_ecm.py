import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.pyplot as plt 

class CellEcm:
    
    def __init__(self, data, params):
        """
        Initialize with HPPC battery cell data and model parameters.
        """
        self.current = data.current
        self.time = data.time
        self.voltage = data.voltage
        self.id1, self.id2, self.id3, self.id4 = data.get_indices()
        
        self.eta_chg = params.eta_chg
        self.eta_dis = params.eta_dis
        self.q_cell = params.q_cell

    
    @staticmethod
    def func_ttc(t, a, b, c, alpha, beta):
        """
        Exponential function for a two time constants model (TTC).
        """
        return a - b * np.exp(-alpha * t) - c * np.exp(-beta * t)

    @staticmethod
    def get_rtau(rctau, z):
        
        # determine index where z is close to soc parameters
        soc = np.arange(0.0, 1.00, 0.05)[::-1]
        idx = abs(soc - z).argmin()

        # return resistor and tau values at z
        tau1 = rctau[:, 0][idx]
        tau2 = rctau[:, 1][idx]
        r0 = rctau[:, 2][idx]
        r1 = rctau[:, 3][idx]
        r2 = rctau[:, 4][idx]
        return tau1, tau2, r0, r1, r2

    def soc(self):
        
        current = self.current
        time = self.time
        q = self.q_cell * 3600
        dt = np.diff(time)

        nc = len(current)
        z = np.ones(nc)

        for k in range(1, nc):
            i = current[k]
            if i < 0:
                eta = self.eta_chg
            else:
                eta = self.eta_dis
            z[k] = z[k - 1] - ((eta * i * dt[k - 1]) / q)

        return z

    def soc2(self, data):
    
        current = data.current
        time = data.time
        q = self.q_cell * 3600
        dt = np.diff(time)

        nc = len(current)
        z = np.ones(nc)

        for k in range(1, nc):
            i = current[k]
            if i < 0:
                eta = self.eta_chg
            else:
                eta = self.eta_dis
            z[k] = z[k - 1] - ((eta * i * dt[k - 1]) / q)
        return z



    def ocv(self, soc, pts=False, vz_pts=None):
               
        nc = len(self.id1)
        n_end=len(self.voltage)
        v_pts=[]
        v_pts = self.voltage[0]
        z_pts=[]
        z_pts = soc[0]
        i_pts=[]
        i_pts = self.current[0]
        t_pts=[]
        t_pts = self.time[0]
        
        v_df = self.voltage
        current_df = self.current
        time_df = self.time
                 
        v_array=np.array(v_df.values)
        current_array=np.array(current_df.values)
        time_array=np.array(time_df.values)
        if pts is True:
            for k in range(1, nc):
                aa=self.id1[k]               
                v_pts = np.append(v_pts, v_array[aa])
                z_pts = np.append(z_pts, soc[aa])
                i_pts = np.append(i_pts, current_array[aa])
                t_pts = np.append(t_pts, time_array[aa])
        
            v_pts = np.append(v_pts, v_array[-1])
            z_pts = np.append(z_pts, soc[-1])
            i_pts = np.append(i_pts, current_array[-1])
            t_pts = np.append(t_pts, time_array[-1])
            ocv = np.interp(soc, z_pts[::-1], v_pts[::-1])
            return ocv, i_pts, t_pts, v_pts, z_pts
            
        elif vz_pts is not None:
            v_pts, z_pts = vz_pts
            ocv = np.interp(soc, z_pts[::-1], v_pts[::-1])
            return ocv
        

    def curve_fit_coeff(self, func, ncoeff):
        
        rest_start= self.id4
        rest_end= self.id1
        nrow = len(self.id4)
        coeff = np.zeros((nrow, ncoeff))

        for i in range(0, nrow):
            start = rest_start[i]
            end = rest_end[i+1]
            t_curve_df = self.time[start:end]
            v_curve_df = self.voltage[start:end]
            t_curve_array=np.array(t_curve_df.values)
            v_curve_array=np.array(v_curve_df.values)
            
            t_scale = t_curve_array[-1] - t_curve_array[0]
            guess = v_curve_array[-1], 0.002, 0.01, 0.001, 0.01
            popt, pcov = curve_fit(func, t_scale, v_curve_array, p0=guess)
            coeff[i] = popt
            
            # _, b, c, alpha, beta = coeff[i]
            # new_time=np.arange(0,3600,1)
            # v_model=(v_curve_array[0]+b+c)-b*np.exp(-alpha*new_time)-c*np.exp(-beta*new_time)
            
            # fig, ax = plt.subplots()
            # ax.plot(t_curve_array-t_curve_array[0],  v_curve_array, marker='.', label='exp')
            # ax.plot(new_time, v_model, label='ecm')            
            # plt.ylim([2.6, 4.5]) 
            # plt.show()
            

        return coeff

    def rctau_ttc(self, coeff):
       
        #id0, id1, id2, _, _, = self.idd
        
        s1=self.id1  #휴식끝
        s2=self.id2  #펄스시작
        s3=self.id3  #펄스끝
        s4=self.id4  #휴식시작
        
        nrow = len(s4)
        rctau = np.zeros((nrow, 7))

        for k in range(nrow):
            di = abs(self.current[s1[k]] - self.current[s2[k]])
            dt = self.time[s1[k+1]] - self.time[s4[k]]
            dv = abs(self.voltage[s1[k]] - self.voltage[s2[k]])

            _, b, c, alpha, beta = coeff[k]

            tau1 = 1 / alpha
            tau2 = 1 / beta
            r0 = dv / di
            r1 = b / ((1 - np.exp(-dt / tau1)) * di)
            r2 = c / ((1 - np.exp(-dt / tau2)) * di)
            c1 = tau1 / r1
            c2 = tau2 / r2

            rctau[k] = tau1, tau2, r0, r1, r2, c1, c2

        return rctau

 


    def vt(self, soc, ocv, rctau):
        """
        Determine voltage from equivalent circuit model.
        """
        dt = np.diff(self.time)     # length of each time step, dt is not constant
        nc = len(self.current)      # total number of time steps based on current
        v0 = np.zeros(nc)           # initialize v0 array
        v1 = np.zeros(nc)           # initialize v1 array
        v2 = np.zeros(nc)           # initialize v2 array

        for k in range(1, nc):
            i = self.current[k]

            # get parameters at state of charge
            tau1, tau2, r0, r1, r2 = self.get_rtau(rctau, soc[k])

            # voltage in r0 resistor
            v0[k] = r0 * i

            # voltage in c1 capacitor
            tm1 = v1[k - 1] * np.exp(-dt[k - 1] / tau1)
            tm2 = r1 * (1 - np.exp(-dt[k - 1] / tau1)) * i
            v1[k] = tm1 + tm2

            # voltage in c2 capacitor
            tm3 = v2[k - 1] * np.exp(-dt[k - 1] / tau2)
            tm4 = r2 * (1 - np.exp(-dt[k - 1] / tau2)) * i
            v2[k] = tm3 + tm4

        vt = ocv - v0 - v1 - v2
        return vt

    

    def vt2(self, data, soc, ocv, rctau):
        """
        Determine voltage from equivalent circuit model.
        """
        dt = np.diff(data.time)     # length of each time step, dt is not constant
        nc = len(data.current)      # total number of time steps based on current
        v0 = np.zeros(nc)           # initialize v0 array
        v1 = np.zeros(nc)           # initialize v1 array
        v2 = np.zeros(nc)           # initialize v2 array
        ocv_new =np.zeros(nc)
        x= soc
        y= ocv
        for k in range(1, nc):
            i = data.current[k]

            # get parameters at state of charge
            tau1, tau2, r0, r1, r2 = self.get_rtau(rctau, soc[k])

            # voltage in r0 resistor
            v0[k] = r0 * i

            # voltage in c1 capacitor
            tm1 = v1[k - 1] * np.exp(-dt[k - 1] / tau1)
            tm2 = r1 * (1 - np.exp(-dt[k - 1] / tau1)) * i
            v1[k] = tm1 + tm2

            # voltage in c2 capacitor
            tm3 = v2[k - 1] * np.exp(-dt[k - 1] / tau2)
            tm4 = r2 * (1 - np.exp(-dt[k - 1] / tau2)) * i
            v2[k] = tm3 + tm4
            
            linear_func=interpolate.interp1d(x, y, kind='linear')
            ocv_new[k]=linear_func(soc[k])
            ocv_new[0]=ocv_new[1]
                        
        vt2 = ocv_new - v0 - v1 - v2
        return vt2

    

   