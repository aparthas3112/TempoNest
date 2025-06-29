#!/usr/bin/env python
"""
Enhanced TempoNest Corner Plot Tool

A simple, flexible tool for creating corner plots from MultiNest chains output.
Automatically detects parameter names from paramnames.txt (if available) or uses
common TempoNest parameter patterns as fallback.

NEW: Comparison Mode - Overlay posteriors from two different analysis results!

Parameter Detection:
    1. Reads paramnames.txt from same directory as chains file (preferred)
    2. Falls back to pattern matching based on number of parameters
    3. Uses generic param_0, param_1, etc. as last resort

Basic Usage:
    # Plot all parameters (auto-detects from paramnames.txt if available)
    python plot_corner.py chains.txt
    
    # Plot specific parameters by index
    python plot_corner.py chains.txt --params 0,1,2
    
    # Plot specific parameters by name (with auto-detection)
    python plot_corner.py chains.txt --params red_amp,red_gamma,efac
    
    # Plot all parameters of a specific type (works with legacy and new formats)
    python plot_corner.py chains.txt --params EFAC    # All EFAC parameters (EFAC1,EFAC2... or efac_group_...)
    python plot_corner.py chains.txt --params EQUAD   # All EQUAD parameters (EQUAD1,EQUAD2... or equad_group_...)
    python plot_corner.py chains.txt --params ECORR   # All ECORR parameters (ECORR1,ECORR2... or ecorr::...)
    
    # Mix individual and group selections
    python plot_corner.py chains.txt --params RedAmp,EFAC,EQUAD
    
    # Include truth values
    python plot_corner.py chains.txt --truth -15.0,3.5,0.1
    
    # Save to file instead of displaying
    python plot_corner.py chains.txt --output corner_plot.png

Comparison Mode:
    # Compare two analyses with same parameter structure
    python plot_corner.py chains1.txt --compare /path/to/dir2/ \
        --param-map-1 red_amp,red_gamma,efac \
        --param-map-2 red_amp,red_gamma,efac
    
    # Compare parameter groups (works with legacy and new formats)
    python plot_corner.py chains1.txt --compare /path/to/dir2/ \
        --param-map-1 ECORR \
        --param-map-2 ECORR
    
    # Compare with different parameter orderings
    python plot_corner.py chains1.txt --compare /path/to/dir2/ \
        --param-map-1 0,1,2 \
        --param-map-2 2,0,1
    
    # Compare with custom labels
    python plot_corner.py chains1.txt --compare /path/to/dir2/ \
        --param-map-1 red_amp,red_gamma,efac \
        --param-map-2 red_amp,red_gamma,efac \
        --label-1 "GPU Analysis" --label-2 "CPU Analysis"
    
    # Multiple comparisons (up to 10 datasets)
    python plot_corner.py main_chains.txt \
        --compare dir2/ --compare dir3/ --compare dir4/ \
        --param-map-1 red_amp,red_gamma,efac \
        --param-map-2 red_amp,red_gamma,efac \
        --param-map-2 red_amp,red_gamma,efac \
        --param-map-2 red_amp,red_gamma,efac \
        --label-1 "Main" --label-2 "Test 1" --label-2 "Test 2" --label-2 "Test 3"
    
    # Different files in comparison directories
    python plot_corner.py post_equal_weights.dat \
        --compare dir2/ --compare dir3/ \
        --compare-file chains.txt --compare-file other_chains.txt \
        --param-map-1 red_amp,efac \
        --param-map-2 0,2 --param-map-2 1,3
        
    # Works with various paramnames file formats:
    # - paramnames.txt, *.paramnames, *paramnames*, paramnames
"""

import argparse
import numpy as np
import sys
import os

# Only import plotting libraries when actually needed
def import_plotting_libraries():
    global plt, corner
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: 'matplotlib' package not found. Install with: pip install matplotlib")
        sys.exit(1)
    
    try:
        import corner
    except ImportError:
        print("Error: 'corner' package not found. Install with: pip install corner")
        sys.exit(1)

# Common TempoNest parameter names and their display labels
TEMPONEST_PARAMS = {
    # Red noise
    'red_amp': r'$\log_{10} A_{red}$',
    'red_amplitude': r'$\log_{10} A_{red}$',
    'redamp': r'$\log_{10} A_{red}$',
    'red_gamma': r'$\gamma_{red}$',
    'red_spectral_index': r'$\gamma_{red}$',
    'redgamma': r'$\gamma_{red}$',
    'redslope': r'$\gamma_{red}$',
    
    # DM noise  
    'dm_amp': r'$\log_{10} A_{DM}$',
    'dm_amplitude': r'$\log_{10} A_{DM}$',
    'dmamp': r'$\log_{10} A_{DM}$',
    'dm_gamma': r'$\gamma_{DM}$',
    'dm_spectral_index': r'$\gamma_{DM}$',
    'dmgamma': r'$\gamma_{DM}$',
    'dmslope': r'$\gamma_{DM}$',
    
    # EFAC/EQUAD
    'efac': r'$\log_{10}$ EFAC',
    'equad': r'$\log_{10}$ EQUAD', 
    'efac_global': r'$\log_{10}$ EFAC',
    'equad_global': r'$\log_{10}$ EQUAD',
    
    # Solar Wind
    'electron_density': r'$n_e^{SW}$ (scaling)',
    'solar_wind_density': r'$n_e^{SW}$ (scaling)',
    'sw_density': r'$n_e^{SW}$ (scaling)',
    'log_amplitude': r'$\log_{10} A_{SW}$',
    'sw_log_amp': r'$\log_{10} A_{SW}$',
    'solar_wind_amp': r'$\log_{10} A_{SW}$',
    'stoch_sw': r'$\log_{10} A_{SW}$',
    
    # Timing parameters
    'f0': r'$\Delta F_0$ (Hz)',
    'f1': r'$\Delta F_1$ (Hz/s)',
    'f2': r'$\Delta F_2$ (Hz/s²)',
    'raj': r'$\Delta$ RA (mas)',
    'decj': r'$\Delta$ DEC (mas)',
    'pmra': r'$\Delta \mu_\alpha$ (mas/yr)',
    'pmdec': r'$\Delta \mu_\delta$ (mas/yr)',
    'px': r'$\Delta$ Parallax (mas)',
    'pb': r'$\Delta P_b$ (days)',
    'a1': r'$\Delta a_1$ (lt-s)',
    'e': r'$\Delta e$',
    'omega': r'$\Delta \omega$ (deg)',
    't0': r'$\Delta T_0$ (MJD)',
    'ecc': r'$\Delta e$',
    'om': r'$\Delta \omega$ (deg)',
    'tasc': r'$\Delta T_{asc}$ (MJD)',
    
    # Other
    'likelihood': 'Log Likelihood',
    'loglike': 'Log Likelihood',
    'weight': 'Weight',
    'prior': 'Prior'
}

def detect_parameter_names(filename, num_params):
    """
    Attempt to detect parameter names from various paramnames files or use defaults.
    """
    chains_dir = os.path.dirname(os.path.abspath(filename))
    
    # Try multiple possible paramnames file patterns
    paramnames_patterns = [
        'paramnames.txt',
        '*.paramnames',
        '*paramnames*',
        'paramnames',
    ]
    
    paramnames_file = None
    for pattern in paramnames_patterns:
        if '*' in pattern:
            # Use glob to find files matching pattern
            import glob
            matches = glob.glob(os.path.join(chains_dir, pattern))
            if matches:
                # Use the first match (could be enhanced to pick best match)
                paramnames_file = matches[0]
                break
        else:
            # Direct file check
            candidate = os.path.join(chains_dir, pattern)
            if os.path.exists(candidate):
                paramnames_file = candidate
                break
    
    if paramnames_file:
        try:
            with open(paramnames_file, 'r') as f:
                param_names = [line.strip() for line in f if line.strip()]
            
            # Verify we have the right number of parameters
            if len(param_names) == num_params:
                print(f"Using parameter names from: {paramnames_file}")
                return param_names
            elif len(param_names) == num_params - 1:
                # Likely missing likelihood column
                print(f"Using parameter names from: {paramnames_file} (adding likelihood)")
                return param_names + ['likelihood']
            else:
                print(f"Warning: {os.path.basename(paramnames_file)} has {len(param_names)} names but chains has {num_params} columns")
                print("Falling back to pattern detection...")
        except Exception as e:
            print(f"Warning: Could not read {paramnames_file}: {e}")
            print("Falling back to pattern detection...")
    
    # Fall back to pattern detection
    # Common TempoNest naming patterns
    common_patterns = [
        # 11-parameter model (red noise + timing + efac)
        ['red_amp', 'red_gamma', 'dm_amp', 'dm_gamma', 'f0', 'f1', 'raj', 'decj', 'pmra', 'pmdec', 'efac'],
        
        # Simple red noise model
        ['red_amp', 'red_gamma', 'efac'],
        
        # Red + DM noise model  
        ['red_amp', 'red_gamma', 'dm_amp', 'dm_gamma', 'efac', 'equad'],
        
        # Solar wind models
        ['electron_density', 'efac'],  # Deterministic solar wind + EFAC
        ['log_amplitude', 'efac'],  # Stochastic solar wind + EFAC
        ['red_amp', 'red_gamma', 'electron_density', 'efac'],  # Red noise + deterministic SW
        ['red_amp', 'red_gamma', 'log_amplitude', 'efac'],  # Red noise + stochastic SW
        ['electron_density', 'log_amplitude', 'efac'],  # Both solar wind components
    ]
    
    # Try to match based on number of parameters (excluding likelihood)
    param_count = num_params - 1 if num_params > 1 else num_params
    
    for pattern in common_patterns:
        if len(pattern) == param_count:
            return pattern + ['likelihood'] if num_params > len(pattern) else pattern
    
    # Default naming scheme
    param_names = [f'param_{i}' for i in range(param_count)]
    if num_params > param_count:
        param_names.append('likelihood')
    
    return param_names

def get_display_labels(param_names):
    """
    Convert parameter names to display labels.
    """
    labels = []
    for name in param_names:
        name_lower = name.lower().strip()
        if name_lower in TEMPONEST_PARAMS:
            labels.append(TEMPONEST_PARAMS[name_lower])
        else:
            # Clean up the name for display
            clean_name = name.replace('_', ' ').title()
            labels.append(clean_name)
    return labels

def parse_parameter_selection(param_str, total_params, param_names):
    """
    Parse parameter selection string (indices, names, or parameter groups).
    
    Special group names:
    - 'EFAC' or 'efac': All EFAC parameters (legacy: EFAC1,EFAC2... or new: efac_group_...)
    - 'EQUAD' or 'equad': All EQUAD parameters (legacy: EQUAD1,EQUAD2... or new: equad_group_...)
    - 'ECORR' or 'ecorr': All ECORR parameters (legacy: ECORR1,ECORR2... or new: ecorr::...)
    """
    if not param_str:
        return list(range(total_params))
    
    selections = [p.strip() for p in param_str.split(',')]
    indices = []
    
    for sel in selections:
        try:
            # Try as index
            idx = int(sel)
            if 0 <= idx < total_params:
                indices.append(idx)
            else:
                print(f"Warning: Index {idx} out of range (0-{total_params-1})")
        except ValueError:
            # Check for parameter group selection
            sel_lower = sel.lower()
            
            if sel_lower in ['efac']:
                # Find all EFAC parameters (both legacy and new formats)
                efac_indices = [i for i, name in enumerate(param_names) 
                              if (name.lower().startswith('efac_') or 
                                  name.lower().startswith('efac') and name.lower()[4:].isdigit())]
                if efac_indices:
                    indices.extend(efac_indices)
                    print(f"Found {len(efac_indices)} EFAC parameters: {[param_names[i] for i in efac_indices]}")
                else:
                    print(f"Warning: No EFAC parameters found")
                continue
                
            elif sel_lower in ['equad']:
                # Find all EQUAD parameters (both legacy and new formats)
                equad_indices = [i for i, name in enumerate(param_names) 
                               if (name.lower().startswith('equad_') or 
                                   name.lower().startswith('equad') and name.lower()[5:].isdigit())]
                if equad_indices:
                    indices.extend(equad_indices)
                    print(f"Found {len(equad_indices)} EQUAD parameters: {[param_names[i] for i in equad_indices]}")
                else:
                    print(f"Warning: No EQUAD parameters found")
                continue
                
            elif sel_lower in ['ecorr']:
                # Find all ECORR parameters (legacy ECORR1/ECORR2, new ecorr::, or old per_backend::)
                ecorr_indices = [i for i, name in enumerate(param_names) 
                               if (name.lower().startswith('ecorr::') or 
                                   name.lower().startswith('per_backend::') or
                                   (name.lower().startswith('ecorr') and name.lower()[5:].isdigit()))]
                if ecorr_indices:
                    indices.extend(ecorr_indices)
                    print(f"Found {len(ecorr_indices)} ECORR parameters: {[param_names[i] for i in ecorr_indices]}")
                else:
                    print(f"Warning: No ECORR parameters found")
                continue
            
            # Try as individual parameter name
            found = False
            for i, name in enumerate(param_names):
                if name.lower() == sel_lower:
                    indices.append(i)
                    found = True
                    break
            if not found:
                print(f"Warning: Parameter '{sel}' not found")
    
    # Remove duplicates while preserving user-specified order
    seen = set()
    ordered_indices = []
    for idx in indices:
        if idx not in seen:
            seen.add(idx)
            ordered_indices.append(idx)
    return ordered_indices

def apply_transformations(data, param_names):
    """
    Apply common transformations (e.g., EFAC from log to linear).
    """
    transformed_data = data.copy()
    transformed_names = param_names.copy()
    
    for i, name in enumerate(param_names):
        name_lower = name.lower()
        # Convert EFAC from log10 to linear scale
        if 'efac' in name_lower and not name_lower.endswith('_linear'):
            transformed_data[:, i] = 10**(data[:, i])
            # Update label
            if name_lower in TEMPONEST_PARAMS:
                transformed_names[i] = 'EFAC'  # Remove log10 from label
    
    return transformed_data, transformed_names

def load_comparison_data(compare_dir, compare_file, input_filename):
    """
    Load comparison data from a second directory.
    """
    if compare_file:
        compare_path = os.path.join(compare_dir, compare_file)
    else:
        # Use same filename as input
        input_basename = os.path.basename(input_filename)
        compare_path = os.path.join(compare_dir, input_basename)
    
    if not os.path.exists(compare_path):
        raise FileNotFoundError(f"Comparison file not found: {compare_path}")
    
    try:
        data = np.loadtxt(compare_path)
        print(f"Loaded comparison data: {data.shape[0]} samples with {data.shape[1]} parameters from {compare_path}")
        
        # Ensure 2D
        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        # Detect parameter names for comparison data
        param_names = detect_parameter_names(compare_path, data.shape[1])
        if len(param_names) != data.shape[1]:
            param_names = [f'param_{i}' for i in range(data.shape[1])]
        
        return data, param_names
        
    except Exception as e:
        raise RuntimeError(f"Error loading comparison file '{compare_path}': {e}")

def map_parameters_for_comparison(param_map_str, total_params, param_names, dataset_name):
    """
    Parse parameter mapping for comparison datasets.
    Supports indices, individual parameter names, and parameter groups (EFAC, EQUAD, ECORR).
    """
    if not param_map_str:
        return list(range(total_params))
    
    selections = [p.strip() for p in param_map_str.split(',')]
    indices = []
    
    for sel in selections:
        try:
            # Try as index
            idx = int(sel)
            if 0 <= idx < total_params:
                indices.append(idx)
            else:
                print(f"Warning: Index {idx} out of range for {dataset_name} (0-{total_params-1})")
        except ValueError:
            # Check for parameter group selection
            sel_lower = sel.lower()
            
            if sel_lower in ['efac']:
                # Find all EFAC parameters (both legacy and new formats)
                efac_indices = [i for i, name in enumerate(param_names) 
                              if (name.lower().startswith('efac_') or 
                                  name.lower().startswith('efac') and name.lower()[4:].isdigit())]
                if efac_indices:
                    indices.extend(efac_indices)
                    print(f"Found {len(efac_indices)} EFAC parameters in {dataset_name}: {[param_names[i] for i in efac_indices]}")
                else:
                    print(f"Warning: No EFAC parameters found in {dataset_name}")
                continue
                
            elif sel_lower in ['equad']:
                # Find all EQUAD parameters (both legacy and new formats)
                equad_indices = [i for i, name in enumerate(param_names) 
                               if (name.lower().startswith('equad_') or 
                                   name.lower().startswith('equad') and name.lower()[5:].isdigit())]
                if equad_indices:
                    indices.extend(equad_indices)
                    print(f"Found {len(equad_indices)} EQUAD parameters in {dataset_name}: {[param_names[i] for i in equad_indices]}")
                else:
                    print(f"Warning: No EQUAD parameters found in {dataset_name}")
                continue
                
            elif sel_lower in ['ecorr']:
                # Find all ECORR parameters (legacy ECORR1/ECORR2, new ecorr::, or old per_backend::)
                ecorr_indices = [i for i, name in enumerate(param_names) 
                               if (name.lower().startswith('ecorr::') or 
                                   name.lower().startswith('per_backend::') or
                                   (name.lower().startswith('ecorr') and name.lower()[5:].isdigit()))]
                if ecorr_indices:
                    indices.extend(ecorr_indices)
                    print(f"Found {len(ecorr_indices)} ECORR parameters in {dataset_name}: {[param_names[i] for i in ecorr_indices]}")
                else:
                    print(f"Warning: No ECORR parameters found in {dataset_name}")
                continue
            
            # Try as individual parameter name
            found = False
            for i, name in enumerate(param_names):
                if name.lower() == sel_lower:
                    indices.append(i)
                    found = True
                    break
            if not found:
                print(f"Warning: Parameter '{sel}' not found in {dataset_name}")
    
    # Remove duplicates while preserving user-specified order
    seen = set()
    ordered_indices = []
    for idx in indices:
        if idx not in seen:
            seen.add(idx)
            ordered_indices.append(idx)
    return ordered_indices

def main():
    parser = argparse.ArgumentParser(
        description="Create corner plots from TempoNest MultiNest chains",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument("input", 
                       help="Input chains file (e.g., chains.txt, post_equal_weights.dat)")
    
    parser.add_argument("--params", "-p", 
                       help="Parameter selection: indices (0,1,2), names (red_amp,efac), or groups (EFAC,EQUAD,ECORR). "
                            "Groups select all parameters of that type. Mix individual and group selections. "
                            "Default: plot all parameters except likelihood")
    
    parser.add_argument("--truth", "-t", 
                       help="Comma-separated truth values (e.g., -15.0,3.5,0.1)")
    
    parser.add_argument("--output", "-o", 
                       help="Output filename (PNG/PDF). Default: display plot")
    
    parser.add_argument("--exclude-likelihood", action="store_true",
                       help="Exclude likelihood column from plotting (default if >10 params)")
    
    parser.add_argument("--include-likelihood", action="store_true", 
                       help="Force include likelihood column in plotting")
    
    parser.add_argument("--transform", action="store_true", default=True,
                       help="Apply standard transformations (EFAC log→linear) (default)")
    
    parser.add_argument("--no-transform", action="store_true",
                       help="Don't apply transformations, plot raw values")
    
    parser.add_argument("--title", 
                       help="Plot title")
    
    parser.add_argument("--dpi", type=int, default=100,
                       help="DPI for saved plots (default: 100)")
    
    # Comparison arguments - support multiple datasets
    parser.add_argument("--compare", "-c", action="append",
                       help="Path to comparison directory. Can be used multiple times for multiple comparisons.")
    
    parser.add_argument("--compare-file", action="append",
                       help="Specific filename in comparison directory (default: same as input file). Can be used multiple times.")
    
    parser.add_argument("--param-map-1", 
                       help="Parameter names/indices for first dataset (main input file)")
    
    parser.add_argument("--param-map-2", action="append",
                       help="Parameter names/indices for comparison dataset. Can be used multiple times.")
    
    # Dataset labels for legend
    parser.add_argument("--label-1", default="Dataset 1",
                       help="Label for first dataset in legend (default: 'Dataset 1')")
    
    parser.add_argument("--label-2", action="append",
                       help="Label for comparison datasets in legend. Can be used multiple times.")
    
    args = parser.parse_args()
    
    # Import plotting libraries only when needed
    import_plotting_libraries()
    
    # Load data
    if not os.path.exists(args.input):
        print(f"Error: File '{args.input}' not found")
        return 1
    
    try:
        data = np.loadtxt(args.input)
        print(f"Loaded {data.shape[0]} samples with {data.shape[1]} parameters")
    except Exception as e:
        print(f"Error loading file '{args.input}': {e}")
        return 1
    
    # Ensure 2D
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    total_params = data.shape[1]
    
    # Detect parameter names
    param_names = detect_parameter_names(args.input, total_params)
    if len(param_names) != total_params:
        param_names = [f'param_{i}' for i in range(total_params)]
    
    print(f"Detected parameters: {param_names}")
    
    # Handle comparison mode - support multiple datasets
    comparison_datasets = []
    comparison_mode = bool(args.compare)
    
    if comparison_mode:
        num_comparisons = len(args.compare)
        
        # Validate arguments
        if not args.param_map_1:
            print("Error: When using --compare, --param-map-1 must be specified")
            return 1
        
        if not args.param_map_2 or len(args.param_map_2) != num_comparisons:
            print(f"Error: Need exactly {num_comparisons} --param-map-2 arguments for {num_comparisons} comparison datasets")
            print("Example: --compare dir1/ --compare dir2/ --param-map-2 params1 --param-map-2 params2")
            return 1
        
        # Load all comparison datasets
        for i, compare_dir in enumerate(args.compare):
            if not os.path.exists(compare_dir):
                print(f"Error: Comparison directory '{compare_dir}' not found")
                return 1
            
            # Get comparison file for this dataset
            compare_file = None
            if args.compare_file and i < len(args.compare_file):
                compare_file = args.compare_file[i]
            
            try:
                comp_data, comp_param_names = load_comparison_data(
                    compare_dir, compare_file, args.input)
                comparison_datasets.append({
                    'data': comp_data,
                    'param_names': comp_param_names,
                    'param_map': args.param_map_2[i],
                    'directory': compare_dir
                })
                print(f"Comparison dataset {i+1} parameters: {comp_param_names}")
            except (FileNotFoundError, RuntimeError) as e:
                print(f"Error loading comparison dataset {i+1}: {e}")
                return 1
    
    # Determine which parameters to plot
    if comparison_mode:
        # Use parameter mappings for comparison
        plot_indices_1 = map_parameters_for_comparison(args.param_map_1, total_params, param_names, "dataset 1")
        
        # Validate all comparison datasets have matching parameter counts
        for i, comp_dataset in enumerate(comparison_datasets):
            comp_indices = map_parameters_for_comparison(
                comp_dataset['param_map'], 
                comp_dataset['data'].shape[1], 
                comp_dataset['param_names'], 
                f"dataset {i+2}"
            )
            comp_dataset['plot_indices'] = comp_indices
            
            if len(plot_indices_1) != len(comp_indices):
                print(f"Error: Parameter mappings must have same length. Dataset 1: {len(plot_indices_1)}, Dataset {i+2}: {len(comp_indices)}")
                return 1
        
        if not plot_indices_1:
            print("Error: No valid parameters selected for comparison")
            return 1
        
        plot_indices = plot_indices_1  # Use first dataset indices for primary data
        
    elif args.params:
        plot_indices = parse_parameter_selection(args.params, total_params, param_names)
    else:
        # Default: exclude likelihood if many parameters or if requested
        plot_indices = list(range(total_params))
        if (total_params > 10 and 'likelihood' in param_names[-1].lower()) or args.exclude_likelihood:
            plot_indices = plot_indices[:-1]
        if args.include_likelihood and 'likelihood' in param_names[-1].lower():
            plot_indices = list(range(total_params))  # Include all
    
    if not plot_indices:
        print("Error: No parameters selected for plotting")
        return 1
    
    if comparison_mode:
        print(f"Plotting parameters from dataset 1: {[param_names[i] for i in plot_indices_1]}")
        for i, comp_dataset in enumerate(comparison_datasets):
            print(f"Plotting parameters from dataset {i+2}: {[comp_dataset['param_names'][j] for j in comp_dataset['plot_indices']]}")
    else:
        print(f"Plotting parameters: {[param_names[i] for i in plot_indices]}")
    
    # Select data and parameter names
    plot_data = data[:, plot_indices]
    plot_param_names = [param_names[i] for i in plot_indices]
    
    # Handle comparison data
    if comparison_mode:
        # Process all comparison datasets
        for comp_dataset in comparison_datasets:
            comp_plot_data = comp_dataset['data'][:, comp_dataset['plot_indices']]
            comp_plot_param_names = [comp_dataset['param_names'][i] for i in comp_dataset['plot_indices']]
            
            # Apply transformations
            if args.transform and not args.no_transform:
                comp_plot_data, comp_plot_param_names = apply_transformations(comp_plot_data, comp_plot_param_names)
            
            # Store processed data back
            comp_dataset['plot_data'] = comp_plot_data
            comp_dataset['plot_param_names'] = comp_plot_param_names
        
        # Apply transformations to main dataset
        if args.transform and not args.no_transform:
            plot_data, plot_param_names = apply_transformations(plot_data, plot_param_names)
    else:
        # Apply transformations to single dataset
        if args.transform and not args.no_transform:
            plot_data, plot_param_names = apply_transformations(plot_data, plot_param_names)
    
    # Get display labels (use labels from first dataset for comparison mode)
    display_labels = get_display_labels(plot_param_names)
    
    # Parse truth values
    truths = None
    if args.truth:
        try:
            truth_values = [float(val.strip()) for val in args.truth.split(',')]
            if len(truth_values) == len(plot_indices):
                truths = truth_values
            else:
                print(f"Warning: {len(truth_values)} truth values provided for {len(plot_indices)} parameters")
                # Pad with None or truncate
                if len(truth_values) < len(plot_indices):
                    truths = truth_values + [None] * (len(plot_indices) - len(truth_values))
                else:
                    truths = truth_values[:len(plot_indices)]
        except ValueError as e:
            print(f"Error parsing truth values: {e}")
    
    # Create corner plot
    if comparison_mode:
        # Define colors for multiple datasets
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        num_datasets = len(comparison_datasets) + 1  # +1 for main dataset
        
        if num_datasets > len(colors):
            print(f"Warning: Too many datasets ({num_datasets}). Only first {len(colors)} will have distinct colors.")
        
        # Create figure with first dataset
        figure = corner.corner(
            plot_data, 
            labels=display_labels,
            truths=truths,
            show_titles=True,
            title_fmt=".3f",
            quantiles=[0.16, 0.5, 0.84],
            smooth=True,
            hist_kwargs={'density': True, 'alpha': 0.7},
            contour_kwargs={'colors': [colors[0]], 'alpha': 0.7},
            color=colors[0]
        )
        
        # Overplot all comparison datasets
        for i, comp_dataset in enumerate(comparison_datasets):
            color_idx = (i + 1) % len(colors)
            corner.corner(
                comp_dataset['plot_data'],
                fig=figure,
                labels=display_labels,
                smooth=True,
                hist_kwargs={'density': True, 'alpha': 0.7},
                contour_kwargs={'colors': [colors[color_idx]], 'alpha': 0.7},
                color=colors[color_idx]
            )
        
        # Create legend with custom labels
        from matplotlib.lines import Line2D
        legend_elements = []
        
        # Main dataset label
        main_label = args.label_1 if args.label_1 != "Dataset 1" else args.label_1
        legend_elements.append(Line2D([0], [0], color=colors[0], alpha=0.7, label=main_label))
        
        # Comparison dataset labels
        for i, comp_dataset in enumerate(comparison_datasets):
            color_idx = (i + 1) % len(colors)
            if args.label_2 and i < len(args.label_2):
                label = args.label_2[i]
            else:
                label = f"Dataset {i+2}"
            legend_elements.append(Line2D([0], [0], color=colors[color_idx], alpha=0.7, label=label))
        
        figure.legend(handles=legend_elements, loc='upper right')
        
    else:
        # Single dataset plot
        figure = corner.corner(
            plot_data, 
            labels=display_labels,
            truths=truths,
            show_titles=True,
            title_fmt=".3f",
            quantiles=[0.16, 0.5, 0.84],
            smooth=True,
            hist_kwargs={'density': True}
        )
    
    # Add title if provided
    if args.title:
        figure.suptitle(args.title, fontsize=16)
    
    # Save or display
    if args.output:
        plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    else:
        plt.show()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())