import streamlit as st
import json
import io
import zipfile
import os

# Initialize update counter if not present.
if "update_counter" not in st.session_state:
    st.session_state.update_counter = 0

# Set the webpage tab title
st.set_page_config(page_title="TempoNest JSON Configurator")

st.title("TempoNest JSON Configurator")
st.markdown("Generate JSON configuration files for TempoNest. Choose to download as a zip file or save directly to a local directory.")

# --- File Naming Settings (Sidebar) ---
st.sidebar.header("File Naming Settings")
main_filename = st.sidebar.text_input("Main Config File Name", value="main.json")
settings_filename = st.sidebar.text_input("Settings File Name", value="settings.json")
sampler_filename = st.sidebar.text_input("Sampler File Name", value="sampler.json")
noise_filename = st.sidebar.text_input("Noise Model File Name", value="noise.json")
timing_filename = st.sidebar.text_input("Timing Model File Name", value="timing.json")
zip_filename = st.sidebar.text_input("Zip File Name", value="tnest_configs.zip")

# --- Output Options (Sidebar) ---
st.sidebar.header("Output Options")
output_format = st.sidebar.selectbox("Output Format", ["Save to Local Directory", "Download Zip"], index=0)

if output_format == "Save to Local Directory":
    output_directory = st.sidebar.text_input("Output Directory Path", 
                                            value="./",
                                            help="Local path where JSON files will be saved")
    folder_name = st.sidebar.text_input("Folder Name", 
                                      value="tnest_config",
                                      help="Name of the folder to create (will be created inside the output directory path)")

# --- Tabs for Configuration Sections ---
tabs = st.tabs(["Main Config", "Settings", "Sampler", "Noise Model", "Timing Model"])

# --- Main Config Tab ---
with tabs[0]:
    st.header("Main Config")
    st.markdown("This configuration will include the names of the other JSON configuration files.")
    
    # Show current output settings
    if output_format == "Download Zip":
        st.info(f"📦 **Output**: Download zip file named `{zip_filename}`")
    else:
        if output_format == "Save to Local Directory":
            output_path = os.path.join(output_directory, folder_name)
            st.info(f"📁 **Output**: Save to directory `{output_path}`")
        else:
            st.info("📁 **Output**: Not configured")

# --- Settings Tab ---
with tabs[1]:
    st.header("Settings")
    st.markdown("Define global parameters for TempoNest.")
    use_original_errors = st.checkbox("Use Original Errors", value=True)
    num_tempo2_its = st.number_input("Number of Tempo2 Iterations", value=1, min_value=1, step=1)
    test_mode = st.checkbox("Test Mode", value=False)

# --- Sampler Tab ---
with tabs[2]:
    st.header("Sampler")
    st.markdown("Configure the sampler settings.")
    sampler_type = st.selectbox("Sampler Type", ["multinest", "polychord"], index=0)
    output_root = st.text_input("Output Root", value="results/TNest-")
    sample = st.checkbox("Sample", value=True)
    importance_sampling = st.selectbox("Importance Sampling", options=[0, 1], index=0)
    constant_efficiency = st.selectbox("Constant Efficiency", options=[0, 1], index=1)
    efficiency = st.number_input("Efficiency", value=0.1, min_value=0.0, max_value=1.0, step=0.01)
    live_points = st.number_input("Live Points", value=500, min_value=1, step=1)

# --- Callback functions for Noise Model ---
def delete_noise_element(index):
    st.session_state.noise_elements_list.pop(index)

def add_noise_element():
    new_noise = {
        "element_name": "Power Law Red Noise",
        "model_type": "Global",  # Default for EFAC/EQUAD
        "days_per_coeff": 30.0,  # Default for Power Law noise models
        "parameters": [
            {"name": "amplitude", "description": "log amplitude of the power law noise process",
             "prior_type": "log_uniform", "include": True, "fit": True, "min_value": -18, "max_value": -10},
            {"name": "spectral_index", "description": "spectral index (red)",
             "prior_type": "uniform", "include": True, "fit": True, "min_value": 0, "max_value": 7}
        ]
    }
    st.session_state.noise_elements_list.append(new_noise)

# --- Noise Model Tab ---
with tabs[3]:
    st.header("Noise Model")
    st.markdown("Define noise model elements. You can delete any element or add a new one.")
    
    st.info("""
    **Noise Model Options**:
    • **Power Law Red/DM Noise**: Stochastic processes with time-span dependent frequency grids
      - `days_per_coeff`: Controls frequency spacing (default: 30 days/coeff)
      - Number of frequencies = floor(observation_span / days_per_coeff)
    • **EFAC/EQUAD**: Error scaling (EFAC) and quadrature addition (EQUAD) - Global or Per-Flag
    • **ECORR**: Epoch correlations for timing noise - Multiple flags with shared priors
      - Epoch window: 10 seconds (legacy compatible)
      - Automatically detects backends for each selected flag
    • **Deterministic Solar Wind**: Corrects systematic electron density variations  
    • **Stochastic Solar Wind**: Accounts for timing noise from solar wind fluctuations
    
    💡 **Solar Wind Requirements**: Tempo2 must provide tdis2 values and NE_SW parameter in .par file.
    📊 **Frequency Grid**: Automatically calculated from observation time span for scientific consistency.
    """)
    
    # Initialize default noise elements if not already in session state.
    if "noise_elements_list" not in st.session_state:
        st.session_state.noise_elements_list = [
            {
                "element_name": "Power Law Red Noise",
                "days_per_coeff": 30.0,
                "parameters": [
                    {"name": "amplitude", "description": "log amplitude of the power law noise process",
                     "prior_type": "log_uniform", "include": True, "fit": True, "min_value": -18, "max_value": -10},
                    {"name": "spectral_index", "description": "spectral index (red)",
                     "prior_type": "uniform", "include": True, "fit": True, "min_value": 0, "max_value": 7}
                ]
            },
            {
                "element_name": "Power Law DM Noise",
                "days_per_coeff": 30.0,
                "parameters": [
                    {"name": "amplitude", "description": "log amplitude of the DM noise process",
                     "prior_type": "log_uniform", "include": True, "fit": True, "min_value": -18, "max_value": -10},
                    {"name": "spectral_index", "description": "spectral index (DM)",
                     "prior_type": "uniform", "include": True, "fit": True, "min_value": 0, "max_value": 7}
                ]
            },
            {
                "element_name": "EFAC",
                "model_type": "Per-Flag",
                "flag": "-group",
                "parameters": [
                    {"name": "per_flag", "description": "EFAC per flag value",
                     "prior_type": "uniform", "include": True, "fit": True, "min_value": -1, "max_value": 0.7, "flag": "-group"}
                ]
            },
            {
                "element_name": "EQUAD",
                "model_type": "Per-Flag", 
                "flag": "-group",
                "parameters": [
                    {"name": "per_flag", "description": "EQUAD per flag value",
                     "prior_type": "log_uniform", "include": True, "fit": True, "min_value": -9, "max_value": -3, "flag": "-group"}
                ]
            }
            # Uncomment below to include solar wind models by default
            # {
            #     "element_name": "Deterministic Solar Wind",
            #     "parameters": [
            #         {"name": "electron_density", "description": "Solar wind electron density scaling factor",
            #          "prior_type": "uniform", "include": True, "fit": True, "min_value": 0.0, "max_value": 10.0}
            #     ]
            # },
            # {
            #     "element_name": "Stochastic Solar Wind", 
            #     "parameters": [
            #         {"name": "log_amplitude", "description": "Log10 amplitude of stochastic solar wind noise",
            #          "prior_type": "log_uniform", "include": True, "fit": True, "min_value": -18.0, "max_value": -10.0}
            #     ]
            # }
        ]

    # Render each noise element in an expander without the element header text.
    for i, elem in enumerate(st.session_state.noise_elements_list):
        with st.expander("", expanded=True):
            old_elem_name = elem["element_name"]
            new_elem_name = st.selectbox("Element Name",
                                         ["Power Law Red Noise", "Power Law DM Noise", "EFAC", "EQUAD", "ECORR", "Deterministic Solar Wind", "Stochastic Solar Wind"],
                                         index=["Power Law Red Noise", "Power Law DM Noise", "EFAC", "EQUAD", "ECORR", "Deterministic Solar Wind", "Stochastic Solar Wind"].index(elem["element_name"]),
                                         key=f"noise_elem_name_{i}")
            st.session_state.noise_elements_list[i]["element_name"] = new_elem_name
            
            # Check if element type changed - if so, reset to default parameters
            element_type_changed = (old_elem_name != new_elem_name)
            
            if new_elem_name in ["Deterministic Solar Wind"]:
                # Deterministic Solar Wind configuration
                st.markdown("**Deterministic Solar Wind**: Corrects timing residuals for systematic variations in solar wind electron density.")
                col1, col2 = st.columns(2)
                with col1:
                    # Use correct defaults for deterministic solar wind (0.0 to 10.0)
                    if element_type_changed or not elem.get("parameters"):
                        current_min = 0.0  # Default for deterministic solar wind
                    else:
                        current_min = elem["parameters"][0].get("min_value", 0.0)
                    min_value = st.number_input("Electron Density Min Value",
                                                value=current_min,
                                                help="Minimum value for electron density scaling factor (typical range: 0.0-10.0)",
                                                key=f"det_sw_min_{i}")
                with col2:
                    # Use correct defaults for deterministic solar wind (0.0 to 10.0)
                    if element_type_changed or not elem.get("parameters"):
                        current_max = 10.0  # Default for deterministic solar wind
                    else:
                        current_max = elem["parameters"][0].get("max_value", 10.0)
                    max_value = st.number_input("Electron Density Max Value",
                                                value=current_max,
                                                help="Maximum value for electron density scaling factor (typical range: 0.0-10.0)",
                                                key=f"det_sw_max_{i}")
                st.session_state.noise_elements_list[i]["parameters"] = [{
                    "name": "electron_density",
                    "description": "Solar wind electron density scaling factor",
                    "prior_type": "uniform",
                    "include": True,
                    "fit": True,
                    "min_value": min_value,
                    "max_value": max_value
                }]
            elif new_elem_name in ["Stochastic Solar Wind"]:
                # Stochastic Solar Wind configuration
                st.markdown("**Stochastic Solar Wind**: Adds frequency-dependent white noise from unmodeled solar wind fluctuations.")
                col1, col2 = st.columns(2)
                with col1:
                    min_value = st.number_input("Log Amplitude Min Value",
                                                value=elem["parameters"][0].get("min_value", -18.0) if elem["parameters"] else -18.0,
                                                help="Minimum value for log10 amplitude prior (typical range: -18.0 to -10.0)",
                                                key=f"stoch_sw_min_{i}")
                with col2:
                    max_value = st.number_input("Log Amplitude Max Value",
                                                value=elem["parameters"][0].get("max_value", -10.0) if elem["parameters"] else -10.0,
                                                help="Maximum value for log10 amplitude prior (typical range: -18.0 to -10.0)",
                                                key=f"stoch_sw_max_{i}")
                st.session_state.noise_elements_list[i]["parameters"] = [{
                    "name": "log_amplitude",
                    "description": "Log10 amplitude of stochastic solar wind noise",
                    "prior_type": "log_uniform",
                    "include": True,
                    "fit": True,
                    "min_value": min_value,
                    "max_value": max_value
                }]
            elif new_elem_name in ["EFAC", "EQUAD"]:
                # Model type selection (Global vs Per-Flag)
                model_type = st.selectbox("Model Type", 
                                        ["Global", "Per-Flag"], 
                                        index=0 if elem.get("model_type", "Per-Flag") == "Global" else 1,
                                        key=f"noise_model_type_{i}")
                st.session_state.noise_elements_list[i]["model_type"] = model_type
                
                if model_type == "Global":
                    col1, col2 = st.columns(2)
                    with col1:
                        min_value = st.number_input("Global Min Value",
                                                    value=elem["parameters"][0].get("min_value", -1 if new_elem_name=="EFAC" else -9),
                                                    key=f"noise_global_min_{i}")
                    with col2:
                        max_value = st.number_input("Global Max Value",
                                                    value=elem["parameters"][0].get("max_value", 0.7 if new_elem_name=="EFAC" else -3),
                                                    key=f"noise_global_max_{i}")
                    st.session_state.noise_elements_list[i]["parameters"] = [{
                        "name": "global",
                        "description": f"global scaling for {new_elem_name}" if new_elem_name=="EFAC" else f"global quadrature term for {new_elem_name}",
                        "prior_type": "uniform" if new_elem_name=="EFAC" else "log_uniform",
                        "include": True,
                        "fit": True,
                        "min_value": min_value,
                        "max_value": max_value
                    }]
                else:  # Per-Flag
                    # Flag selection
                    flag_name = st.text_input("Flag Name", 
                                            value=elem.get("flag", "-group"), 
                                            help="Tempo2 flag to use (e.g., -group, -fe, -sys, -f)",
                                            key=f"noise_flag_{i}")
                    st.session_state.noise_elements_list[i]["flag"] = flag_name
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        min_value = st.number_input("Per-Flag Min Value",
                                                    value=elem["parameters"][0].get("min_value", -1 if new_elem_name=="EFAC" else -9),
                                                    key=f"noise_perflag_min_{i}")
                    with col2:
                        max_value = st.number_input("Per-Flag Max Value",
                                                    value=elem["parameters"][0].get("max_value", 0.7 if new_elem_name=="EFAC" else -3),
                                                    key=f"noise_perflag_max_{i}")
                    st.session_state.noise_elements_list[i]["parameters"] = [{
                        "name": "per_flag",
                        "description": f"{new_elem_name} per flag value",
                        "prior_type": "uniform" if new_elem_name=="EFAC" else "log_uniform",
                        "include": True,
                        "fit": True,
                        "min_value": min_value,
                        "max_value": max_value,
                        "flag": flag_name
                    }]
            elif new_elem_name == "ECORR":
                # ECORR configuration with multiple flags and shared priors
                st.markdown("**ECORR**: Epoch correlations with multiple flag support. Select multiple flags and set shared priors.")
                
                # Flag selection (multiple flags)
                available_flags = ["-B", "-fe", "-f", "-sys", "-chan", "-group"]
                current_flags = elem.get("flags", ["-B"])  # Default to -B flag
                
                selected_flags = st.multiselect(
                    "Select ECORR Flags", 
                    options=available_flags,
                    default=current_flags,
                    help="Select one or more Tempo2 flags. Each flag creates separate ECORR elements with shared priors.",
                    key=f"ecorr_flags_{i}"
                )
                st.session_state.noise_elements_list[i]["flags"] = selected_flags
                
                if not selected_flags:
                    st.warning("⚠️ Please select at least one flag for ECORR")
                else:
                    st.info(f"📊 Will create {len(selected_flags)} ECORR element(s) for flags: {', '.join(selected_flags)}")
                
                # Shared prior configuration
                st.markdown("**Shared Priors** (applied to all selected flags):")
                col1, col2 = st.columns(2)
                with col1:
                    min_value = st.number_input("ECORR Min Value (log10)",
                                              value=elem.get("min_value", -9),
                                              help="Minimum log10(ECORR) value in microseconds",
                                              key=f"ecorr_min_{i}")
                with col2:
                    max_value = st.number_input("ECORR Max Value (log10)",
                                              value=elem.get("max_value", -3),
                                              help="Maximum log10(ECORR) value in microseconds", 
                                              key=f"ecorr_max_{i}")
                
                # Advanced configuration
                st.markdown("**Advanced ECORR Settings:**")
                col3, col4 = st.columns(2)
                with col3:
                    epoch_window = st.number_input("Epoch Window (seconds)",
                                                 value=elem.get("epoch_window", 10.0),
                                                 min_value=0.1,
                                                 help="Time window for grouping TOAs into epochs (default: 10 seconds, legacy compatible)",
                                                 key=f"ecorr_window_{i}")
                with col4:
                    min_toas = st.number_input("Min TOAs per Epoch",
                                             value=elem.get("min_toas_per_epoch", 1),
                                             min_value=1,
                                             help="Minimum number of TOAs required to form an epoch (default: 1, legacy compatible)",
                                             key=f"ecorr_min_toas_{i}")
                
                # Store configuration
                st.session_state.noise_elements_list[i]["min_value"] = min_value
                st.session_state.noise_elements_list[i]["max_value"] = max_value
                st.session_state.noise_elements_list[i]["epoch_window"] = epoch_window
                st.session_state.noise_elements_list[i]["min_toas_per_epoch"] = min_toas
                
                # Set parameter structure (will be expanded during JSON generation)
                st.session_state.noise_elements_list[i]["parameters"] = [{
                    "name": "per_backend",
                    "description": "ECORR per backend (shared priors across flags)",
                    "prior_type": "log_uniform",
                    "include": True,
                    "fit": True,
                    "min_value": min_value,
                    "max_value": max_value
                }]
            else:
                # Add days_per_coeff configuration
                days_per_coeff = st.number_input("Days per Coefficient", 
                                                value=elem.get("days_per_coeff", 30.0), 
                                                min_value=1.0, 
                                                step=1.0,
                                                help="Time span per frequency coefficient in days (default: 30.0). Number of frequencies = floor(time_span/days_per_coeff)",
                                                key=f"noise_days_per_coeff_{i}")
                st.session_state.noise_elements_list[i]["days_per_coeff"] = days_per_coeff
                
                col1, col2 = st.columns(2)
                with col1:
                    amp_include = st.checkbox("Include amplitude", value=elem["parameters"][0].get("include", True), key=f"noise_amp_include_{i}")
                    amp_fit = st.checkbox("Fit amplitude", value=elem["parameters"][0].get("fit", True), key=f"noise_amp_fit_{i}")
                with col2:
                    amp_min = st.number_input("Amplitude Min Value", value=elem["parameters"][0].get("min_value", -18), key=f"noise_amp_min_{i}")
                    amp_max = st.number_input("Amplitude Max Value", value=elem["parameters"][0].get("max_value", -10), key=f"noise_amp_max_{i}")
                col3, col4 = st.columns(2)
                with col3:
                    si_include = st.checkbox("Include spectral_index", value=elem["parameters"][1].get("include", True), key=f"noise_si_include_{i}")
                    si_fit = st.checkbox("Fit spectral_index", value=elem["parameters"][1].get("fit", True), key=f"noise_si_fit_{i}")
                with col4:
                    si_min = st.number_input("spectral_index Min Value", value=elem["parameters"][1].get("min_value", 0), key=f"noise_si_min_{i}")
                    si_max = st.number_input("spectral_index Max Value", value=elem["parameters"][1].get("max_value", 7), key=f"noise_si_max_{i}")
                st.session_state.noise_elements_list[i]["parameters"] = [
                    {"name": "amplitude", "description": "log amplitude of the power law noise process", "prior_type": "log_uniform",
                     "include": amp_include, "fit": amp_fit, "min_value": amp_min, "max_value": amp_max},
                    {"name": "spectral_index", "description": "spectral index (red)", "prior_type": "uniform",
                     "include": si_include, "fit": si_fit, "min_value": si_min, "max_value": si_max}
                ]
            st.button("Delete this Noise Element", key=f"noise_delete_{i}", on_click=delete_noise_element, args=(i,))
    st.button("Add Noise Element", on_click=add_noise_element)

    noise_elements_for_json = []
    for elem in st.session_state.noise_elements_list:
        if elem["element_name"] == "ECORR":
            # Handle ECORR specially - create separate elements for each flag
            selected_flags = elem.get("flags", [])
            for flag in selected_flags:
                json_elem = {
                    "name": "ECORR",
                    "parameters": [{
                        "name": "per_backend",
                        "prior_type": "log_uniform",
                        "include": True,
                        "fit": True,
                        "min_value": elem.get("min_value", -9),
                        "max_value": elem.get("max_value", -3),
                        "flag": flag,
                        "epoch_window": elem.get("epoch_window", 10.0) / 86400.0,  # Convert to days
                        "min_toas_per_epoch": elem.get("min_toas_per_epoch", 1)
                    }]
                }
                noise_elements_for_json.append(json_elem)
        else:
            # Handle other elements normally
            json_elem = {
                "name": elem["element_name"],
                "parameters": elem["parameters"]
            }
            # Add days_per_coeff for Power Law noise models
            if elem["element_name"] in ["Power Law Red Noise", "Power Law DM Noise"] and "days_per_coeff" in elem:
                json_elem["days_per_coeff"] = elem["days_per_coeff"]
            noise_elements_for_json.append(json_elem)


# --- Callback functions for Timing Model ---
def delete_timing_param(index):
    st.session_state.timing_params.pop(index)

def add_timing_param():
    new_param = {
        "name": f"PARAM{len(st.session_state.timing_params)+1}",
        "prior_type": "uniform",
        "include": True,
        "fit": True
    }
    st.session_state.timing_params.append(new_param)

# --- Timing Model Tab ---
with tabs[4]:
    st.header("Timing Model")
    st.markdown("Define timing model parameters. All parameters will be grouped under a single 'Timing Model' element.")
    
    # Timing Model Mode Selection
    st.subheader("Timing Model Configuration")
    
    timing_mode = st.radio(
        "Select timing model handling mode:",
        options=[
            "Marginalise over all timing parameters", 
            "Fit for all timing model parameters", 
            "Manual timing parameter selection"
        ],
        index=0,  # Default to marginalisation
        help="Choose how TempoNest should handle timing model parameters from your .par file"
    )
    
    if timing_mode == "Marginalise over all timing parameters":
        st.info("🔄 **Marginalisation Mode**: All fitted parameters from your .par file will be marginalised analytically. This is the fastest and most robust option for most analyses.")
        marginalise_all = True
        fit_all = False
        sigma_multiplier = 10  # Not used but set for consistency
        
    elif timing_mode == "Fit for all timing model parameters":
        st.info("📊 **Fit All Mode**: All fitted parameters from your .par file will be included in sampling with configurable uncertainty ranges. This provides full parameter posteriors but is computationally expensive.")
        marginalise_all = False
        fit_all = True
        
        sigma_multiplier = st.number_input(
            "Sigma multiplier for parameter ranges",
            value=10.0,
            min_value=1.0,
            max_value=50.0,
            step=1.0,
            help="Prior ranges will be: fitted_value ± (sigma_multiplier × fitted_uncertainty)"
        )
        st.caption(f"Parameter ranges will be ±{sigma_multiplier}σ around fitted values")
        
    else:  # Manual timing parameter selection
        st.info("⚙️ **Manual Mode**: You can explicitly select which timing parameters to sample below. All other fitted parameters from your .par file will be automatically marginalised.")
        marginalise_all = False
        fit_all = False
        sigma_multiplier = 10  # Not used but set for consistency
    
    # Initialize timing parameters if not already in session state.
    if "timing_params" not in st.session_state:
        st.session_state.timing_params = [
            {"name": "RAJ", "prior_type": "uniform", "include": False, "fit": True},
            {"name": "DECJ", "prior_type": "uniform", "include": False, "fit": True},
            {"name": "F0", "prior_type": "uniform", "include": False, "fit": True},
            {"name": "F1", "prior_type": "uniform", "include": False, "fit": True}
        ]
    
    # Show parameter configuration based on selected mode
    if timing_mode == "Marginalise over all timing parameters":
        st.subheader("Parameter Overrides (Advanced)")
        expander_context = st.expander("Click to override marginalisation for specific parameters", expanded=False)
        show_params = True
    elif timing_mode == "Fit for all timing model parameters":
        st.subheader("Automatic Parameter Configuration")
        st.markdown("All fitted parameters from your .par file will be automatically included with the specified sigma multiplier. No manual configuration needed.")
        show_params = False
        expander_context = None
    else:  # Manual timing parameter selection
        st.subheader("Manual Parameter Configuration")
        st.markdown("Configure specific timing parameters to be included in sampling. All other fitted parameters will be marginalised.")
        show_params = True
        expander_context = None
    
    # Context manager for the parameters section
    if show_params and expander_context is not None:
        with expander_context:
            st.markdown("Add specific timing parameters to override marginalisation and include them in sampling instead.")
            # Render each timing parameter
            for i, param in enumerate(st.session_state.timing_params):
                container = st.container()
                with container:
                    col1, col2, col3, col4, col5 = st.columns([2, 2, 1, 1, 1])
                    new_name = col1.text_input("Parameter Name", value=param["name"], key=f"timing_name_{i}")
                    new_prior = col2.selectbox("Prior Type", ["uniform", "log_uniform"],
                                                index=0 if param["prior_type"]=="uniform" else 1,
                                                key=f"timing_prior_{i}")
                    new_include = col3.checkbox("Include", value=param["include"], key=f"timing_include_{i}")
                    new_fit = col4.checkbox("Fit", value=param["fit"], key=f"timing_fit_{i}")
                    # Use on_click callback to delete immediately.
                    col5.button("Delete", key=f"timing_delete_{i}", on_click=delete_timing_param, args=(i,))
                    
                    # Show uncertainty multiplier inputs only if parameter is included
                    if new_include:
                        st.markdown(f"**{new_name} Prior Range (uncertainty multipliers):**")
                        col_min, col_max = st.columns(2)
                        with col_min:
                            new_min = st.number_input("Min Multiplier", value=param.get("min_value", -5.0), 
                                                    help="Prior min = fitted_value + (this × fitted_uncertainty)", 
                                                    key=f"timing_min_{i}")
                        with col_max:
                            new_max = st.number_input("Max Multiplier", value=param.get("max_value", 5.0), 
                                                    help="Prior max = fitted_value + (this × fitted_uncertainty)", 
                                                    key=f"timing_max_{i}")
                        st.markdown("---")
                    else:
                        new_min = param.get("min_value", -5.0)
                        new_max = param.get("max_value", 5.0)
                    
                    # Update session state with any edits.
                    st.session_state.timing_params[i] = {
                        "name": new_name,
                        "prior_type": new_prior,
                        "include": new_include,
                        "fit": new_fit,
                        "min_value": new_min,
                        "max_value": new_max
                    }
            st.button("Add Timing Parameter", on_click=add_timing_param)
    elif show_params and expander_context is None:
        # Render each timing parameter (not in expander)
        for i, param in enumerate(st.session_state.timing_params):
            container = st.container()
            with container:
                col1, col2, col3, col4, col5 = st.columns([2, 2, 1, 1, 1])
                new_name = col1.text_input("Parameter Name", value=param["name"], key=f"timing_name_{i}")
                new_prior = col2.selectbox("Prior Type", ["uniform", "log_uniform"],
                                            index=0 if param["prior_type"]=="uniform" else 1,
                                            key=f"timing_prior_{i}")
                new_include = col3.checkbox("Include", value=param["include"], key=f"timing_include_{i}")
                new_fit = col4.checkbox("Fit", value=param["fit"], key=f"timing_fit_{i}")
                # Use on_click callback to delete immediately.
                col5.button("Delete", key=f"timing_delete_{i}", on_click=delete_timing_param, args=(i,))
                
                # Show uncertainty multiplier inputs only if parameter is included
                if new_include:
                    st.markdown(f"**{new_name} Prior Range (uncertainty multipliers):**")
                    col_min, col_max = st.columns(2)
                    with col_min:
                        new_min = st.number_input("Min Multiplier", value=param.get("min_value", -5.0), 
                                                help="Prior min = fitted_value + (this × fitted_uncertainty)", 
                                                key=f"timing_min_{i}")
                    with col_max:
                        new_max = st.number_input("Max Multiplier", value=param.get("max_value", 5.0), 
                                                help="Prior max = fitted_value + (this × fitted_uncertainty)", 
                                                key=f"timing_max_{i}")
                    st.markdown("---")
                else:
                    new_min = param.get("min_value", -5.0)
                    new_max = param.get("max_value", 5.0)
                
                # Update session state with any edits.
                st.session_state.timing_params[i] = {
                    "name": new_name,
                    "prior_type": new_prior,
                    "include": new_include,
                    "fit": new_fit,
                    "min_value": new_min,
                    "max_value": new_max
                }
        st.button("Add Timing Parameter", on_click=add_timing_param)

    # Group all parameters into one Timing Model element.
    timing_element = {"name": "Timing Model"}
    
    # Configure based on selected timing mode
    if timing_mode == "Marginalise over all timing parameters":
        timing_element["marginalise"] = "all"
        # Only include parameters that have include=True (parameter overrides)
        included_params = []
        for param in st.session_state.timing_params:
            if param["include"]:
                included_params.append({
                    "name": param["name"],
                    "prior_type": param["prior_type"],
                    "include": param["include"],
                    "fit": param["fit"],
                    "min_value": param["min_value"],
                    "max_value": param["max_value"]
                })
        
        # For pure marginalisation mode, only include actual parameter overrides
        timing_element["parameters"] = included_params
        
    elif timing_mode == "Fit for all timing model parameters":
        timing_element["marginalise"] = False
        timing_element["fit_all"] = True
        timing_element["sigma_multiplier"] = sigma_multiplier
        # Provide dummy parameter template for C++ implementation
        timing_element["parameters"] = [{
            "name": "F0",
            "prior_type": "uniform", 
            "include": False,
            "fit": True,
            "min_value": -sigma_multiplier,
            "max_value": sigma_multiplier
        }]
        
    else:  # Manual timing parameter selection
        timing_element["marginalise"] = "manual"
        # Include all parameters with include=True
        included_params = []
        for param in st.session_state.timing_params:
            if param["include"]:
                included_params.append({
                    "name": param["name"],
                    "prior_type": param["prior_type"],
                    "include": param["include"],
                    "fit": param["fit"],
                    "min_value": param["min_value"],
                    "max_value": param["max_value"]
                })
        
        # Always include parameters array for manual mode
        timing_element["parameters"] = included_params if included_params else st.session_state.timing_params
    
    timing_elements = [timing_element]


# --- Generate JSON Files (in Main Config tab) ---
with tabs[0]:
    button_text = "Generate and Download Zip" if output_format == "Download Zip" else "Generate and Save to Directory"
    
    if st.button(button_text):
        # Build the timing model configuration here to access the timing mode variables
        timing_element = {"name": "Timing Model"}
        
        # Configure based on selected timing mode
        if timing_mode == "Marginalise over all timing parameters":
            timing_element["marginalise"] = "all"
            # Only include parameters that have include=True (parameter overrides)
            included_params = []
            for param in st.session_state.timing_params:
                if param["include"]:
                    included_params.append({
                        "name": param["name"],
                        "prior_type": param["prior_type"],
                        "include": param["include"],
                        "fit": param["fit"],
                        "min_value": param["min_value"],
                        "max_value": param["max_value"]
                    })
            
            # For pure marginalisation mode, only include actual parameter overrides
            timing_element["parameters"] = included_params
            
        elif timing_mode == "Fit for all timing model parameters":
            timing_element["marginalise"] = False
            timing_element["fit_all"] = True
            timing_element["sigma_multiplier"] = sigma_multiplier
            # Provide dummy parameter template for C++ implementation
            timing_element["parameters"] = [{
                "name": "F0",
                "prior_type": "uniform", 
                "include": False,
                "fit": True,
                "min_value": -sigma_multiplier,
                "max_value": sigma_multiplier
            }]
            
        else:  # Manual timing parameter selection
            timing_element["marginalise"] = "manual"
            # Include all parameters with include=True
            included_params = []
            for param in st.session_state.timing_params:
                if param["include"]:
                    included_params.append({
                        "name": param["name"],
                        "prior_type": param["prior_type"],
                        "include": param["include"],
                        "fit": param["fit"],
                        "min_value": param["min_value"],
                        "max_value": param["max_value"]
                    })
            
            # Always include parameters array for manual mode
            timing_element["parameters"] = included_params if included_params else st.session_state.timing_params
        
        timing_elements_final = [timing_element]
        
        # Build the configuration dictionaries
        main_config = {
            "INCLUDE": [settings_filename, sampler_filename, noise_filename, timing_filename]
        }
        settings_config = {
            "globals": {
                "use_original_errors": use_original_errors,
                "num_tempo2_its": num_tempo2_its,
                "test_mode": test_mode
            }
        }
        sampler_config = {
            "sampler": {
                "type": sampler_type,
                "output_root": output_root,
                "sample": sample,
                "importance_sampling": importance_sampling,
                "constant_efficiency": constant_efficiency,
                "efficiency": efficiency,
                "live_points": live_points
            }
        }
        noise_model_config = {
            "elements": noise_elements_for_json
        }
        timing_model_config = {
            "elements": timing_elements_final
        }
        
        # Convert each config to a formatted JSON string
        main_json = json.dumps(main_config, indent=4)
        settings_json = json.dumps(settings_config, indent=4)
        sampler_json = json.dumps(sampler_config, indent=4)
        noise_json = json.dumps(noise_model_config, indent=4)
        timing_json = json.dumps(timing_model_config, indent=4)
        
        if output_format == "Download Zip":
            # Create a zip file in memory containing all JSON files
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                zip_file.writestr(main_filename, main_json)
                zip_file.writestr(settings_filename, settings_json)
                zip_file.writestr(sampler_filename, sampler_json)
                zip_file.writestr(noise_filename, noise_json)
                zip_file.writestr(timing_filename, timing_json)
            zip_buffer.seek(0)
            
            st.download_button(
                "Download Zip File",
                data=zip_buffer,
                file_name=zip_filename,
                mime="application/zip"
            )
            st.success("Zip file generated!")
            
        else:  # Save to Local Directory
            try:
                # Create the full path
                full_path = os.path.join(output_directory, folder_name)
                
                # Create directory if it doesn't exist
                os.makedirs(full_path, exist_ok=True)
                
                # Write each JSON file to the directory
                files_written = []
                
                with open(os.path.join(full_path, main_filename), 'w') as f:
                    f.write(main_json)
                    files_written.append(main_filename)
                
                with open(os.path.join(full_path, settings_filename), 'w') as f:
                    f.write(settings_json)
                    files_written.append(settings_filename)
                
                with open(os.path.join(full_path, sampler_filename), 'w') as f:
                    f.write(sampler_json)
                    files_written.append(sampler_filename)
                
                with open(os.path.join(full_path, noise_filename), 'w') as f:
                    f.write(noise_json)
                    files_written.append(noise_filename)
                
                with open(os.path.join(full_path, timing_filename), 'w') as f:
                    f.write(timing_json)
                    files_written.append(timing_filename)
                
                st.success(f"✅ Successfully saved {len(files_written)} JSON files to: `{full_path}`")
                st.info(f"📁 Files created: {', '.join(files_written)}")
                
            except Exception as e:
                st.error(f"❌ Failed to save files: {str(e)}")
                st.info("💡 Make sure the output directory path exists and you have write permissions.")
