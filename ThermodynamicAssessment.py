import numpy as np
import pandas as pd
gamma_C= 1.4 
gamma_H = 1.3
c_PC = 1004.5 #[J/(kg-K)]
c_PH = 1243.66 #[J/kg-K]
R= 287 #[J/kg-K]
fuel_density = 6.843  #[lbs/gallon] #WANT TO CONVERT THIS TO SI
inlet_duct_diam = 3.5 #[inches]

def lbm_to_kg(psia):
    return psia * 0.4536
def lbf_to_N(lbf):
    return 4.4482*lbf
def psi_to_pa(psi):
    return 6894.8 * psi
def celsius_to_kelvin(celsius):
    return celsius + 273.15
def inch_h20_to_pascals(inch_h20):
    return inch_h20*249.1
def rpm_to_rad_per_sec(rpm):
    return (np.pi/30)*rpm
def gal_per_hour_to_kg_per_s(fuel_flow):
    """
    for Jet-A1 fuel
    """
    return (1/3600)*(6.843)*(0.4536) * fuel_flow

def fuel_flow_fraction(m_dot_fuel, m_dot):
    """
    unitless
    """
    return m_dot_fuel/m_dot
def ad_compressor_efficiency(P_t3, P_o, T_t3, T_o):
    """
    Appendix D: Thermodynamic Assesment
    Adiabatic compressor efficiency (unitless)
    Input:
        - P_t3: Compressor Discharge Pressure (Pa)
        - P_o: Cell Pressure [Pa]
        - T_t3: Compressor Discharge Temp [k]
        - T_o: Inlet temperature (T_T1) [k]
    """
    pressure_ratio = P_t3/P_o
    temperature_ratio = T_t3/T_o
    # print((pressure_ratio)**((gamma_C-1)/gamma_C) -1)
    # print(temperature_ratio-1)
    n_c = ((pressure_ratio)**((gamma_C-1)/gamma_C) -1)/(temperature_ratio-1)
    return n_c
def turbine_inlet_temp(T_t5, T_t3, T_o, f):
    """
    Appendix D: Thermodynamic Assesment
    Calculates Turbine Inlet Temperature
    Turbine inlet temperature (T_t4)
    
    """
    c_P_ratio = c_PC/c_PH
    T_t4 = T_t5 + (c_P_ratio)*(T_t3-T_o)/(f+1)
    return T_t4
def isentropic_turbine_exit_temp(T_t5, T_t3, T_o, f, n_c, delta_n):
    """
    Appendix D: Thermodynamic Assesment
    Calculates Isentropic Turbine Exit Temperature [K]
    Parameters:
        -T_t5: Exhaust gas temperature
        -T_t3: Compressor discharge temperature
        -T_o: Inlet temperature (T_t1)
        -f: fuel flow fraction
    n_c: adiabatic compressor efficency;
    delta_n: ?????
    """
    c_P_ratio = c_PC/c_PH
    return T_t5 + c_P_ratio * ((T_t3-T_o)/(f+1))*(1-1/(n_c+delta_n))
def jet_velocity(T_t5, T_t5_s, T_t4, P_o, P_t3):
    """
    Appendix D: Thermodynamic Assesment
    Calculates jet velocity c6:
    Parameters:
        -T_t5_s???? #bruh
        -T_t5: Exhuast gas temperature [K]
        -T_t4: Turbine inlet temperature [K]
        -P_o: cell pressure [Pa]
        -P_t3: Compressor discharge temperature [Pa]
    Return jet velocity
    """
    temp_ratio = T_t4/T_t5_s
    pressure_ratio = P_o/P_t3
    return np.sqrt(np.abs(2*c_PH*T_t5*(1-temp_ratio*(pressure_ratio)**((gamma_H-1)/gamma_H))))
def specific_thrust(T_o, c_6, f):
    """
    Part of Appendix C: Direct Asessment
    Calculates spefic thrust F/(mdot*a_o)
    Params:
        -F: Thrust [N]
    """
    a_o = np.sqrt(gamma_C*R*T_o)
    
    return (f+1)*c_6/a_o

def thrust_specific_fuel_consumption(f, c_6):
    return f/(c_6*(f+1))
def display_info(title, headers, data_lists, spacing_width=20):
    # Determine the number of rows and columns
    num_rows = len(data_lists[0])    
    # Print the title
    print("------------------------------------")
    print(title)
    
    # Print the headers
    header_str = "".join([f"{header:<{spacing_width}}" for header in headers])
    print(header_str)
    
    test_conditions = ["Pre-Test", "Idle Speed", "Speed 1", "Speed 2", "Max Speed", "Post-Test"]
    # Print the data

    for i in range(num_rows):
        condition = test_conditions[i] if i < len(test_conditions) else ""
        row_str = f"{condition:<{spacing_width}}" + "".join([f"{round(data_list[i], 3):<{spacing_width}}" for data_list in data_lists])
        print(row_str)
class ThermodynamicAssessment():
    def __init__(self, df):
        self.compr_effic = []
        self.T_t4 = [] #inlet mass flow
        self.T_t5S = []
        self.c_6 = []
        self.specific_thrusts = []
      
        self.collected_data = {
            "rpms":df["RPMs"].to_numpy(), 
            "cell_pressures": df["Cell Pressures (psia)"].to_numpy(), 
            "inlet_dynamic_pressures": df["Inlet Dynamic Pressures (in H2O)"].to_numpy(), 
            "inlet_temps": df["Inlet Temperatures (°C)"].to_numpy(),
            "compr_discharge_static_pressures": df["Compressor Discharge Static Pressures (psia)"].to_numpy(),
            "compr_discharge_skin_temperatures": df["Compressor Discharge Skin Temperatures (°C)"].to_numpy(),
            "thrusts": df["Thrusts (N)"].to_numpy(),
            "fuel_flows": df["Fuel Flows (kg/s)"].to_numpy(),
            "exhaust_gas_temps": df["Exhaust Gas Temperatures (°C)"].to_numpy(),
            }
        self.rpms = rpm_to_rad_per_sec(self.collected_data["rpms"])
        self.cell_pressures = psi_to_pa(self.collected_data["cell_pressures"])
        self.inlet_dynamic_pressures = inch_h20_to_pascals(self.collected_data["inlet_dynamic_pressures"])
        self.inlet_temps = celsius_to_kelvin(self.collected_data["inlet_temps"])
        self.compr_discharge_static_pressures = psi_to_pa(self.collected_data["compr_discharge_static_pressures"])
        self.compr_discharge_skin_temperatures = celsius_to_kelvin(self.collected_data["compr_discharge_skin_temperatures"])
        self.thrusts = lbf_to_N(self.collected_data["thrusts"])
        self.fuel_flows = gal_per_hour_to_kg_per_s(self.collected_data["fuel_flows"])
        self.exhaust_gas_temps = celsius_to_kelvin(self.collected_data["exhaust_gas_temps"])
        self.inlet_duct_static_pressures = [cell_pressure - inlet_dynamic_pressure for (cell_pressure, inlet_dynamic_pressure) in zip(self.cell_pressures,self.inlet_dynamic_pressures)]#P1
        self.mdots = np.array([0.0048176240735977, 0.052766225144284404, 0.10225058401549976, 0.18673348885890587, 0.21733295816596004, 0.004808650420085111])
    def ad_compressor_efficiency(self, P_t3, P_o, T_t3, T_o):
        """
        Appendix D: Thermodynamic Assesment
        Adiabatic compressor efficiency (unitless)
        Input:
            - P_t3: Compressor Discharge Pressure (Pa)
            - P_o: Cell Pressure [Pa]
            - T_t3: Compressor Discharge Temp [k]
            - T_o: Inlet temperature (T_T1) [k]
        """
        pressure_ratio = P_t3/P_o
        temperature_ratio = T_t3/T_o
        # print((pressure_ratio)**((gamma_C-1)/gamma_C) -1)
        # print(temperature_ratio-1)
        n_c = ((pressure_ratio)**((gamma_C-1)/gamma_C) -1)/(temperature_ratio-1)
        return n_c
    def turbine_inlet_temp(self, T_t5, T_t3, T_o, f):
        """
        Appendix D: Thermodynamic Assesment
        Calculates Turbine Inlet Temperature
        Turbine inlet temperature (T_t4)
        
        """
        c_P_ratio = c_PC/c_PH
        T_t4 = T_t5 + (c_P_ratio)*(T_t3-T_o)/(f+1)
        return T_t4
    def isentropic_turbine_exit_temp(self, T_t5, T_t3, T_o, f, n_c, delta_n):
        """
        Appendix D: Thermodynamic Assesment
        Calculates Isentropic Turbine Exit Temperature [K]
        Parameters:
            -T_t5: Exhaust gas temperature
            -T_t3: Compressor discharge temperature
            -T_o: Inlet temperature (T_t1)
            -f: fuel flow fraction
        n_c: adiabatic compressor efficency;
        delta_n: ?????
        """
        c_P_ratio = c_PC/c_PH
        return T_t5 + c_P_ratio * ((T_t3-T_o)/(f+1))*(1-1/(n_c+delta_n))
    def jet_velocity(self, T_t5, T_t5_s, T_t4, P_o, P_t3):
        """
        Appendix D: Thermodynamic Assesment
        Calculates jet velocity c6:
        Parameters:
            -T_t5_s????
            -T_t5: Exhuast gas temperature [K]
            -T_t4: Turbine inlet temperature [K]
            -P_o: cell pressure [Pa]
            -P_t3: Compressor discharge temperature [Pa]
        Return jet velocity
        """
        temp_ratio = T_t4/T_t5_s
        pressure_ratio = P_o/P_t3
        return np.sqrt(np.abs(2*c_PH*T_t5*(1-temp_ratio*(pressure_ratio)**((gamma_H-1)/gamma_H))))
    def specific_thrust(self, T_o, c_6, f):
        """
        Part of Appendix C: Direct Asessment
        Calculates spefic thrust F/(mdot*a_o)
        Params:
            -F: Thrust [N]
        """
        a_o = np.sqrt(gamma_C*R*T_o)
        
        return (f+1)*c_6/a_o

    def thrust_specific_fuel_consumption(self, f, c_6):
        return f/(c_6*(f+1)) 
            # inlet_diameter = inches_to_meters(inlet_diameter)
            # self.A_1 = calc_inlet_area(inlet_diameter) #INLET AREA
            # self.A_1s = self.A_1*np.ones(6)
    def evaluate(self):
        self.compr_effic = self.ad_compressor_efficiency(self.compr_discharge_static_pressures,self.cell_pressures,self.compr_discharge_skin_temperatures,self.inlet_temps)
        self.compr_effic[0] = 1 #dont want to eval 0/0
        self.compr_effic[5] = 1
        
        self.inlet_temps_SI = np.array(self.inlet_temps)
        self.f = self.fuel_flows/self.mdots
        self.T_t4 = self.turbine_inlet_temp(self.exhaust_gas_temps,self.compr_discharge_skin_temperatures,self.inlet_temps,self.f)
        self.T_t5S = self.isentropic_turbine_exit_temp(self.exhaust_gas_temps,self.compr_discharge_skin_temperatures,self.inlet_temps,self.f,self.compr_effic,0)
        self.c_6 = self.jet_velocity(self.exhaust_gas_temps,self.T_t5S, self.T_t4, self.cell_pressures,self.compr_discharge_static_pressures)
        self.c_6[0] = 0 #Manually set to 0 for Idle, Pre, and Post-Test, since preheating/residual heat messed up the calculation
        self.c_6[1] = 0
        self.c_6[5] = 0

        self.specific_thrusts = self.specific_thrust(self.inlet_temps,self.c_6,self.f)
    def display_info(self, title, headers, data_lists, spacing_width=20):
        # Determine the number of rows and columns
        num_rows = len(data_lists[0])
        print("------------------------------------")
        print(title)
        
        header_str = "".join([f"{header:<{spacing_width}}" for header in headers])
        print(header_str)
        
        test_conditions = ["Pre-Test", "Idle Speed", "Speed 1", "Speed 2", "Max Speed", "Post-Test"]

        for i in range(num_rows):
            condition = test_conditions[i] if i < len(test_conditions) else ""
            row_str = f"{condition:<{spacing_width}}" + "".join([f"{round(data_list[i], 3):<{spacing_width}}" for data_list in data_lists])
            print(row_str)
    def display(self):
        print("FOR THERMODYNAMIC ASSESSMENT:")
        self.display_info("Measured Data:", ["", "P0[Pa]", "T_t1[K]", "P_t3[Pa]", 'T_t3[K]', 'T_t5[K]', "f[-]"], [self.cell_pressures, self.inlet_temps, self.compr_discharge_static_pressures, self.compr_discharge_skin_temperatures, self.exhaust_gas_temps, self.f], spacing_width=20)
        display_info("Estimated Parameters:", ["", "n_c [-]", "T_t4 [K]", "T_t5S [K]", "c_6 [m/s]", "Specific Thrust [-]"], [self.compr_effic, self.T_t4, self.T_t5S, self.c_6,self.specific_thrusts])
if __name__ == "__main__":
    df = pd.read_csv("AssessmentData.csv")
    results = ThermodynamicAssessment(df)
    results.evaluate()
    results.display()
   