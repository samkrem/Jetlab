import numpy as np
import scipy.optimize 

def degrees_to_radians(degrees):
    return degrees * np.pi/180
def radians_to_degrees(radians):
    return radians*180/np.pi
def inches_to_meters(inches):
    return inches * 0.0254
def rpm_to_rad_per_sec(rpm):
    return (np.pi/30)*rpm
def newtons_to_lbf(newtons):
    return newtons/4.4482216152605
def kg_to_lbm(kg):
    return 2.20462 * kg
def sfc_metric_to_imperial(sfc_metric):
    kg_to_lbm = 2.20462
    newtons_to_lbf = 1/4.4482216152605
    sec_to_hrs = 1/3600
    return sfc_metric * kg_to_lbm/(sec_to_hrs*newtons_to_lbf)
class GeometricJetModel():
    """
    Modeling real brayton cycle jet engine for engine operating with normal air
    """
    M_CHOKE = 1 #always one at choke point
    R = 287 #R_T = R_C
    GAMMA_C=1.4 #specifc heat ratio of compressor

    def __init__(self, d_1, d_2, d_gauge, n_blades, beta_1, omega, 
                  T_inf = 288.15, P_inf = 101325, gamma_T = 1.3,
                  EGT = 953.15, c_theta_2=0, n_C = 0.7, 
                  n_T = 0.7, h_f = 43.15*10**6
                  ):
        #Initial Properties of jet engine presented and ideal assumptions, but could be adjusted for jet engine running similar cycle
        self.T_inf = T_inf 
        self.P_inf = P_inf
        self.gamma_T = gamma_T #specific heat ratio

        self.EGT = EGT #exit gas temperature (kelvin)
        self.c_theta_2 = c_theta_2 #theta component of turbine exit absolute velocity(m/s)
        self.n_C = n_C #comp ad eff, educated assumption based on research
        self.n_T = n_T #turb ad eff, educated assumption based on research
        self.h_f = h_f #Fuel heating value (J)
        
        self.d_1 = d_1 #inner diameter of turbine inches
        self.d_2 = d_2 #outer diameter of turbine inches

        self.r_1 = inches_to_meters(d_1) / 2 # m
        self.r_2 = inches_to_meters(d_2) / 2 # m
        self.r_avg = (self.r_1 + self.r_2) / 2 # m
        
        self.d_gauge = inches_to_meters(d_gauge) # m
        self.n_blades = n_blades
        self.beta_1 = beta_1 # degrees
        self.omega = rpm_to_rad_per_sec(omega) # rad/s
    def calc_cp_cv(self, gamma):
        C_p = self.R/(1 - gamma ** (-1))
        C_v = self.R/(gamma - 1)
        return (C_p, C_v)
    def calc_alpha_1(self): #yes
        """
        Calculate alpha_1, stator blade exit angle (radians)
        """
        num = (self.r_2 - self.r_1) * self.n_blades * self.d_gauge
        den = (np.pi * (self.r_2 ** 2 - self.r_1 ** 2))
        self.alpha_1 = np.arccos(num / den)
        return self.alpha_1
    def calc_c_theta_1(self): 
        """       
        Calculates c_theta_1, the theta component of the stator exit absolute velocity (m/s)
        """
        self.c_theta_1 = (self.omega * self.r_avg * np.tan(self.calc_alpha_1())) / (np.tan(self.calc_alpha_1()) - np.tan(self.beta_1))
        return self.c_theta_1

    def calc_work_coefficient(self):
        """
        Calculates work coefficient
        """
        self.work_coefficient = (self.calc_c_theta_1() - self.c_theta_2)/(self.omega * self.r_avg)
        return self.work_coefficient
    def euler_turbine_equation(self): #yes
        """
        Calculates the change in stagnation enthalpy, ht1- ht2 which is just the shaft work of our JetJoe engine
        """
        self.w_shaft = self.omega * (self.r_1 * self.calc_c_theta_1() - self.r_2 * self.c_theta_2)
        return self.w_shaft
    def calc_pressure_ratio(self, T,gamma, n_ad): #yes
        """
        Calculates pressure ratio, for both turbine and compressor its P_high/P_low
        """
        work_spec = self.euler_turbine_equation()
        C_P = self.calc_cp_cv(gamma)[0]       
        return (work_spec * n_ad / (T * C_P) + 1) ** (gamma / (gamma - 1))
    def calc_pi_C(self):
        """
        Calculate pressure ratio for compressor
        """
        self.pi_C = self.calc_pressure_ratio(self.T_inf ,self.GAMMA_C, self.n_C)
        return self.pi_C
    def calc_pi_T(self):
        """
        Uses same pressure ratio equation for compressor except, the efficiency is in denom, and the pressure ratio is inversed
        """
        self.pi_T = 1/self.calc_pressure_ratio(self.EGT,self.gamma_T, 1/self.n_T) #p_high/p_low
        return self.pi_T
    def calc_turbine_inlet_stag_temp(self): #yes
        """
        Calculates turbine inlet stagnation temperature (T_t4)
        """
        C_P_T = self.calc_cp_cv(self.gamma_T)[0]
        self.T_t4 = self.calc_c_theta_1() * self.omega * self.r_avg / C_P_T + self.EGT
        return self.T_t4 

    def calc_fuel_frac(self): #yes
        """
            Calculates fuel fraction
            T_t3 = (1+f)*C_P_T*(T_t4-T_t5)/C_P_C+T_t2
            T_t3 = -(f*h_f/(1+f)+T_t4*C_P_T)/C_P_C
            solve for f which makes both eqs equal (we dont need T_t3 we just put T_t3 in terms of F)
        """
        C_P_T, C_P_C = self.calc_cp_cv(self.gamma_T)[0], self.calc_cp_cv(self.GAMMA_C)[0]
        T_t2,T_t5, T_t4 = self.T_inf, self.EGT, self.calc_turbine_inlet_stag_temp()       
        f_eq = lambda f: (-(1+f)*C_P_T*(T_t4 - T_t5) / C_P_C - T_t2) - ((f * self.h_f / (1 - f) - T_t4 * C_P_T) / C_P_C)
        self.f = scipy.optimize.fsolve(f_eq,0.1)[0] #fuel fraction
        return self.f
    def calc_P_t4(self):
        self.P_t4 = self.P_inf * self.calc_pi_C()
        return self.P_t4
    def calc_NGV_area(self): 
        """
        Calculate area of nozzle guide vain
        """
        self.A_4 = (self.r_2 - self.r_1)*(self.n_blades) * self.d_gauge
        return self.A_4
    def calc_mdot_T(self): #yes
        """
        Calculate mass flow of turbine
        """
        A_4 = self.calc_NGV_area()
        
        P_t4 = self.calc_P_t4()
        T_t4 = self.calc_turbine_inlet_stag_temp()
     
        self.mdot_T =  self.M_CHOKE * A_4 * P_t4 * np.sqrt(self.gamma_T) / (np.sqrt(self.R * T_t4) * ((1 + (self.gamma_T - 1)/2) * self.M_CHOKE**2)**(1 / 2 * (self.gamma_T + 1) / (self.gamma_T - 1)))
        return self.mdot_T
    def calc_mdot_C(self):
        """
        Calculate mass flow of compressor
        """
        self.mdot_C = self.calc_mdot_T()/(1 + self.calc_fuel_frac())
        return self.mdot_C
    def calc_P_t5(self):
        """
        Calculate stagnation pressure at station 5
        """
       
        self.P_t5 =  self.calc_pi_T() * self.calc_P_t4()
        return  self.calc_pi_T() * self.calc_P_t4()
    def calc_M_exit(self):
        """
        Calculates, m_6, exit mach number
        """
        P_t5 = self.calc_P_t5()
        self.M_6 = np.sqrt(2 / (self.gamma_T - 1) * ((P_t5/self.P_inf) ** ((self.gamma_T - 1) / self.gamma_T) - 1))
        
        return self.M_6
    def calc_T_t6(self):
        """
        Calculate stagnation temperature at station 6
        """
        return self.EGT
    def calc_T_6(self):
        T_t6 = self.calc_T_t6()
        M_6 = self.calc_M_exit()
        self.T_6 = T_t6 / (1 + (self.gamma_T - 1) * (M_6 ** 2) / 2)
        return self.T_6
    def calc_a(self, gamma, T):
        """
        Parameters:
            -T: Temperature
            -R: ideal gas constant
            -gamma: specific heat ratio of gas
        
        """
        return np.sqrt(gamma * self.R * T)
    def calc_a_6(self):
        """
        Calculates the speed of sound at station 6
        """
        self.a_6 = self.calc_a(self.gamma_T, self.calc_T_6())
        return self.a_6
    def calc_c_6(self):
        """
        Calcualtes velocity of fluid at station 6
        """
        self.c_6 = self.calc_M_exit() * self.calc_a_6()  
        return self.c_6
    def calc_F(self):
        """
        Calculates the force generated by the engine
        """
        return self.calc_c_6() * self.calc_mdot_T()
    def calc_SFC(self, option = "metric"):
        """
        Calculates specific fuel consumption

        """
        # print("f",self.calc_fuel_frac())
        self.SFC = self.calc_fuel_frac()/(self.calc_c_6() * (1+self.calc_fuel_frac()))
        if option != "metric":
            return sfc_metric_to_imperial(self.SFC)
        else:
            return self.SFC
    def calc_cycle_pr(self):
        """
        Calculates cycle pressure ratio
        """
        self.cycle_PR = self.calc_P_t4()/self.P_inf
        return self.cycle_PR
    def calc_ST(self):
        """
        Calculates the specific thrust
        """

        
        F = self.calc_F()
        a_0 = self.calc_a(self.T_inf, self.GAMMA_C)
        # #np.sqrt(GAMMA_C*R_C*T_inf)
        self.ST = F/(self.calc_mdot_C() * a_0)
        return self.ST
    def model(self):
        """
        Function: Given these input parameters in constructor, calculates information for a JetJoe 1400 Engine
        """
        return ([
                ("Inner Diameter of Turbine:             " ,self.d_1, "       inches"), 
                ("Outer Diameter of Turbine:             ", self.d_2, "       inches"),
                ("Depth of Stator Blade:                 ", self.d_gauge,        "    inches"), 
                ("Number of Stator Blades:               ", self.n_blades,     "         unitless"),
                ("Turbine Rotor Inlet Angle:             ", self.beta_1,     "         degrees"),
                ("Shaft rotation rate:                   ", self.omega,            "   RPM")

                ],

                [
                ("Cycle Pressure Ratio:                  ", self.calc_cycle_pr(),         "      unitless"), 
                ("Inlet Massflow:                        ", self.calc_mdot_C(),           "     kg/s"),
                ("Turbine Inlet Stagnation Temperature:  ", self.calc_turbine_inlet_stag_temp(),             "       K"), 
                ("Work Coefficient:                      ", self.calc_work_coefficient(), "     unitless"),
                ("Specific Thrust:                       ", self.calc_ST(),              "      unitless"), #make dimensionless
                ("Specific Fuel Consumption:             ", self.calc_SFC(),       "  (kg/s)/N"),
                ("Specific Fuel Consumption:             ", self.calc_SFC(option = "imperial"),     "      (lbm/hr)/lbf")
                ]
            
            )
#self.cycle_pr, self.mdot_C, self.calc_turbine_inlet_stag_temp, self.ST, self.SFC, 

    # return (C_p, C_v)
    def display_results(self):
        inputs, results = self.model()
        # inputs, results = values
        print("Inputs")
        print("-----------------------------------------")
        print("Measurement                             Value         Units")
        for input in inputs:
            print(f"{input[0]} {input[1]:.3g} {input[2]}")
        print("\nOutputs")
        print("-----------------------------------------")
        print("Quantity                                Value       Units     ")
        for result in results:
            print(f"{result[0]} {result[1]:.4g} {result[2]}")
if __name__ == "__main__":
    #rules: input in inches, degrees, output calculates in radians, and metric but display is in degrees and meters
    d1_turbine = 1.442 #inches
    d2_turbine = 2.165 #inches
    #area of turbine is ring
    d_exit_nozzle = 1.641 #inches
    d_gauge = 0.225 #inches 

    n_blades = 13 #unitless
    beta_1 = 31 #degrees
    omega = 160000 #RPM
    M_choke = 1 #unitless, occurs in NGV, assume subsonic
    EGT = 953.15 #K, 
    P_inf = 100000 #Pa
    T_inf = 288.15 #K

    jet_engine = GeometricJetModel(d1_turbine, d2_turbine, d_gauge, n_blades, beta_1, omega)
    res = jet_engine.model()
    jet_engine.display_results()
   