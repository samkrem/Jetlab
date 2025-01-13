import numpy as np
import pandas as pd
gamma_C= 1.4 
gamma_H = 1.3
c_PC = 1004.5 #[J/(kg-K)]
c_PH = 1243.66 #[J/kg-K]
R= 287 #[J/kg-K]
fuel_density = 6.843  #[lbs/gallon] #WANT TO CONVERT THIS TO SI
inlet_duct_diam = 3.5 #[inches]'
#Conversion factors
def inches_to_meters(inches):
    return inches * 0.0254
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
    for fuel density
    """
    return gal_per_hour_to_lb_per_hour(fuel_flow)/(2.205*3600) 
def gal_per_hour_to_lb_per_hour(fuel_flow):
    return fuel_density*fuel_flow
def calc_inlet_area(inlet_diameter):
    return np.pi * (inlet_diameter/2)**2

class DirectAssessment():
    def __init__(self, df, inlet_diameter):
        

        self.M1s = []
        self.mdots = [] #inlet mass flow
        self.specific_thrusts = []
        self.SFCs = []
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

        inlet_diameter = inches_to_meters(inlet_diameter)
        self.A_1 = calc_inlet_area(inlet_diameter) #INLET AREA
        self.A_1s = self.A_1*np.ones(6)
    def inlet_mach_no(self,P_o, P_1):
        """
        Part of Appendix C: Direct Assesment
        Caculates inlet mach number M1 [-]
        Params:
            - P_o: cell pressure [Pa]
            - P_1: Inlet duct wall pressure [Pa]
        """
        pressure_ratio = P_o/P_1
        return np.sqrt(
            (2/(gamma_C-1))*(pressure_ratio**((gamma_C-1)/gamma_C)-1)
        )
    def inlet_mass_flow(self, A_1, P_o, T_o, M_1):  
        """
        Part of Appendix C: Direct Asessment
        Calculates inlet mass flow mdot [kg/s]
        Params: 
            - A_1: Inlet area [m**2]
            - P_o: Cell pressure [Pa]
            - T_o: inlet tempreature [K]
            - M_1: inlet mach number [-]
        """
        exp = ((gamma_C+1)/(2*(gamma_C-1)))
        denom = (1 + (M_1**2 * (gamma_C-1))/2)**exp
    
        return (A_1*P_o)/(np.sqrt(T_o)) * np.sqrt(gamma_C/R)*M_1/denom
    def specific_thrust(self,F, mdot, T_o):
        """
        Part of Appendix C: Direct Asessment
        Calculates spefic thrust F/(mdot*a_o)
        Params:
            -F: Thrust [N]
        """
        a_o = np.sqrt(gamma_C*R*T_o)
        return F/(mdot*a_o)

    def specific_fuel_consumption(self,mdot_f, F):
        """
        Part of Appendix C: Direct Asessment
        Calculates specific fuel consumption SFC [(lbm/hr)/lbf]
        [(lbm/hr)/lbf]
        Input: 
            - mdot_f: mass flow of fuel [kg/s]
            - F: thrust
        """
        return mdot_f/F
    def max_PR(self):
        return self.collected_data["compr_discharge_static_pressures"][4]/ self.collected_data["cell_pressures"][4]
    def evaluate(self):
        self.M1s = self.inlet_mach_no(self.cell_pressures, self.inlet_duct_static_pressures)
        self.mdots = self.inlet_mass_flow(self.A_1, self.cell_pressures,self.inlet_temps,self.M1s)
        self.specific_thrusts = self.specific_thrust(self.thrusts,self.mdots, self.inlet_temps)
        self.SFCs = self.specific_fuel_consumption(gal_per_hour_to_lb_per_hour(self.collected_data["fuel_flows"]), self.collected_data["thrusts"])
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
        print("FOR DIRECT ASSESSMENT:")
        self.display_info("Collected Data:", ["", "RPMs[rpm]", "P0[psia]", 'P0-P1[in-h20]', 'T0[°C]', 'P_t3[psia]', 'T_t3[°C])', 'F[lbf]', 'mdot_f[gal/hr]', 'EGT[°C]'], [self.collected_data["rpms"], self.collected_data["cell_pressures"], self.collected_data["inlet_dynamic_pressures"], self.collected_data["inlet_temps"], self.collected_data["compr_discharge_static_pressures"], self.collected_data["compr_discharge_skin_temperatures"], self.collected_data["thrusts"], self.collected_data["fuel_flows"], self.collected_data["exhaust_gas_temps"]], spacing_width=20)
        self.display_info("Measurement Data", ["", "Cell Presure [Pa]", "Inlet Temperatures[K]", "Thrust[N]", "Inlet Area [m^2]", "Inlet Pressure [Pa]"], [self.cell_pressures, self.inlet_temps, self.thrusts,self.A_1s, self.inlet_duct_static_pressures], spacing_width = 28)
        self.display_info("Estimated Params", ["", "M1[-]", "Mass flow [kg/s]", "Specific Thrust [-]", "SFC [(lbf/hr)/lbf]"], [self.M1s, self.mdots, self.specific_thrusts, self.SFCs], spacing_width = 28)

# Sample data
if __name__ == "__main__":
    df = pd.read_csv("AssessmentData.csv")
    results = DirectAssessment(df, inlet_duct_diam)
    results.evaluate()
    results.display()

