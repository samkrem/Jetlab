import pandas as pd
from GeometricJetModel import GeometricJetModel
from ThermodynamicAssessment import ThermodynamicAssessment
from DirectAssessment import DirectAssessment

# Function to calculate percent error between two values
def percent_error(estimated, actual):
    return abs(100 * (estimated - actual) / actual)

# Function to display collected data in a formatted table
def display_info(title, headers, data_lists, spacing_width=20):
    print("-" * 40)
    print(title)

    # Print headers with appropriate spacing
    print("".join(f"{header:<{spacing_width}}" for header in headers))

    # Define test condition labels for rows
    test_conditions = ["Pre-Test", "Idle Speed", "Speed 1", "Speed 2", "Max Speed", "Post-Test"]
    num_rows = len(data_lists[0])  # Determine the number of rows

    # Print each row of data
    for i in range(num_rows):
        condition = test_conditions[i] if i < len(test_conditions) else ""  # Use test condition if available
        row = [f"{condition:<{spacing_width}}"] + [
            f"{round(data[i], 3):<{spacing_width}}" for data in data_lists
        ]
        print("".join(row))

# Function to display results in a formatted table
def display_results(data_lists, spacing_width=25):
    print("-" * 40)

    # Define column units for display
    units = [
        "", "Press. Ratio (-)", "Inlet Mass Flow (kg/s)", 
        "Turbine Inlet Temp. (K)", "Spec. Thrust (-)", 
        "Spec. Fuel Consump. (lbm/hr)/lbf"
    ]
    print("".join(f"{unit:<{spacing_width}}" for unit in units))

    # Define variable names for rows
    variables = [
        "Model", "Thermo Assessment", "Direct Assessment",
        "% Diff w/ Model & Direct", "% Diff w/ Model & Thermo"
    ]

    # Transpose data for easier row-wise display
    transposed_data = list(zip(*data_lists))

    # Print each row of results
    for i, row in enumerate(transposed_data):
        row_str = f"{variables[i]:<{spacing_width}}" + "".join(
            f"{round(value, 3):<{spacing_width}}" if isinstance(value, (int, float)) else f"{value:<{spacing_width}}"
            for value in row
        )
        print(row_str)

if __name__ == "__main__":
    # Input parameters for jet engine model
    d1_turbine, d2_turbine = 1.442, 2.165  # Turbine diameters (inches)
    d_exit_nozzle, d_gauge = 1.641, 0.225  # Nozzle and gauge diameters (inches)
    inlet_duct_diam = 3.5  # Inlet duct diameter (inches)
    n_blades, beta_1, omega = 13, 31, 160000  # Blade count, blade angle (degrees), rotational speed (RPM)
    M_choke, EGT, P_inf, T_inf = 1, 953.15, 100000, 288.15  # Choke Mach, EGT (K), atmospheric pressure (Pa), temperature (K)

    # Read data from CSV file
    df = pd.read_csv("AssessmentData.csv")
    collected_data = {
        "rpms": df["RPMs"].to_numpy(),
        "cell_pressures": df["Cell Pressures (psia)"].to_numpy(),
        "inlet_dynamic_pressures": df["Inlet Dynamic Pressures (in H2O)"].to_numpy(),
        "inlet_temps": df["Inlet Temperatures (°C)"].to_numpy(),
        "compr_discharge_static_pressures": df["Compressor Discharge Static Pressures (psia)"].to_numpy(),
        "compr_discharge_skin_temperatures": df["Compressor Discharge Skin Temperatures (°C)"].to_numpy(),
        "thrusts": df["Thrusts (N)"].to_numpy(),
        "fuel_flows": df["Fuel Flows (kg/s)"].to_numpy(),
        "exhaust_gas_temps": df["Exhaust Gas Temperatures (°C)"].to_numpy(),
    }

    # Display collected data
    display_info(
        "Collected Data:",
        ["", "RPMs[rpm]", "P0[psia]", "P0-P1[in-h20]", "T0[°C]", "P_t3[psia]", "T_t3[°C]", "F[lbf]", "mdot_f[gal/hr]", "EGT[°C]"],
        [collected_data[k] for k in collected_data],
        spacing_width=20
    )

    # Perform thermodynamic and direct assessments
    thermo = ThermodynamicAssessment(df)  # Thermodynamic assessment
    thermo_results = thermo.evaluate()
    thermo.display()

    direct = DirectAssessment(df, inlet_duct_diam)  # Direct assessment
    direct_results = direct.evaluate()
    direct.display()

    # Perform jet engine modeling
    jet_engine = GeometricJetModel(d1_turbine, d2_turbine, d_gauge, n_blades, beta_1, omega)
    predictions = jet_engine.model()
    jet_engine.display_results()

    # Prepare results for display
    PRs = [
        jet_engine.cycle_PR, "-", direct.max_PR(),
        percent_error(jet_engine.cycle_PR, direct.max_PR()), "-"
    ]
    inlet_mass_flows = [
        jet_engine.mdot_C, "-", direct.mdots[4],
        percent_error(jet_engine.mdot_C, direct.mdots[4]), "-"
    ]
    inlet_temps = [
        jet_engine.T_t4, thermo.T_t4[4], "-", "-",
        percent_error(jet_engine.T_t4, thermo.T_t4[4])
    ]
    STs = [
        jet_engine.ST, thermo.specific_thrusts[4], direct.specific_thrusts[4],
        percent_error(jet_engine.ST, direct.specific_thrusts[4]),
        percent_error(jet_engine.ST, thermo.specific_thrusts[4])
    ]
    SFCs = [
        jet_engine.calc_SFC("imperial"), "-", direct.SFCs[4],
        percent_error(jet_engine.calc_SFC("imperial"), direct.SFCs[4]), "-"
    ]

    # Display results
    display_results([PRs, inlet_mass_flows, inlet_temps, STs, SFCs])