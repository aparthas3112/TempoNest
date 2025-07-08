#!/bin/bash
# TempoNest Build Script
# Comprehensive build script for the reorganized TempoNest

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to initialize build environment (autonomous operation)
initialize_build_environment() {
    print_status "Initializing autonomous build environment..."
    
    # Clear potentially problematic environment variables
    # This ensures the build script works regardless of user's .zshrc/.bashrc
    unset MULTINEST_DIR POLYCHORD_DIR ARRAYFIRE_PATH
    unset LDFLAGS  # Reset to clean state
    
    # Set clean working directory paths (will be set properly by setup functions)
    export MULTINEST_DIR=""
    export POLYCHORD_DIR=""
    export ARRAYFIRE_PATH=""
    export LDFLAGS=""
    
    # Preserve essential paths
    export PATH="$PATH"
    
    print_success "Build environment initialized (clean slate)"
}

# Function to auto-detect or download dependencies
setup_dependencies() {
    print_status "Setting up dependencies..."
    
    # Auto-detect Tempo2
    detect_tempo2
    
    # Setup MultiNest (download if needed)
    setup_multinest
    
    # Setup PolyChord (download if needed)
    setup_polychord
    
    # Setup ArrayFire (download if needed, optional)
    setup_arrayfire
}

# Function to detect Tempo2
detect_tempo2() {
    print_status "Detecting Tempo2..."
    
    # Try environment variable first
    if [ -n "$TEMPO2" ] && [ -f "$TEMPO2/bin/tempo2" ]; then
        print_success "TEMPO2 found via environment: $TEMPO2"
        return 0
    fi
    
    # Try command in PATH
    if command_exists tempo2; then
        local tempo2_bin=$(which tempo2)
        export TEMPO2=$(dirname $(dirname $tempo2_bin))
        print_success "TEMPO2 auto-detected: $TEMPO2"
        return 0
    fi
    
    # Try common installation paths
    for path in /usr/local/tempo2 /opt/tempo2 /usr/local /opt/pulsar/tempo2; do
        if [ -f "$path/bin/tempo2" ]; then
            export TEMPO2="$path"
            print_success "TEMPO2 found at: $TEMPO2"
            return 0
        fi
    done
    
    print_error "Tempo2 not found!"
    print_status "Please install Tempo2:"
    print_status "  - Ubuntu: sudo apt-get install tempo2"
    print_status "  - Module: module load tempo2"
    print_status "  - Manual: export TEMPO2=/path/to/tempo2"
    exit 1
}

# Function to setup MultiNest (bundled with repo)
setup_multinest() {
    print_status "Setting up bundled MultiNest..."
    
    local multinest_dir="$(pwd)/external/MultiNest"
    
    # Check if MultiNest source directory exists
    if [ ! -d "$multinest_dir" ]; then
        print_error "MultiNest source directory not found at: $multinest_dir"
        print_status "This suggests the git repository is incomplete"
        exit 1
    fi
    
    # Check if source files exist
    if [ ! -f "$multinest_dir/nested.F90" ]; then
        print_error "MultiNest source files not found in: $multinest_dir"
        print_status "Expected to find nested.F90 and other source files"
        exit 1
    fi
    
    # Set the directory - build_multinest() will handle compilation
    export MULTINEST_DIR="$multinest_dir"
    
    # Check if already compiled
    if [ -f "$multinest_dir/libnest3.so" ] || [ -f "$multinest_dir/libnest3.a" ]; then
        print_success "MultiNest source found with existing libraries (libnest3)"
    elif [ -f "$multinest_dir/lib/libmultinest.a" ] || [ -f "$multinest_dir/lib/libmultinest.so" ]; then
        print_success "MultiNest source found with existing libraries (libmultinest)"
    else
        print_success "MultiNest source found - will compile during build"
    fi
}

# Function to setup ArrayFire (optional)
setup_arrayfire() {
    print_status "Setting up ArrayFire (optional for GPU acceleration)..."
    
    # Check if ArrayFire already exists in various possible locations
    local possible_existing_dirs=(
        "$(pwd)/external/ArrayFire/ArrayFire-3.9.0-Linux"
        "$(pwd)/external/ArrayFire-3.9.0-Linux"
        "$(pwd)/external/ArrayFire"
    )
    
    for arrayfire_dir in "${possible_existing_dirs[@]}"; do
        if [ -d "$arrayfire_dir" ] && ([ -f "$arrayfire_dir/lib64/libaf.so" ] || [ -f "$arrayfire_dir/lib/libaf.so" ]); then
            print_success "ArrayFire already installed at: $arrayfire_dir"
            # Ensure lib directory exists to prevent libtool errors
            if [ ! -d "$arrayfire_dir/lib" ]; then
                mkdir -p "$arrayfire_dir/lib"
                print_status "Created missing lib directory in ArrayFire"
            fi
            export ARRAYFIRE_PATH="$arrayfire_dir"
            return 0
        fi
    done
    
    # Look for ArrayFire installer in scripts directory
    local arrayfire_installer="$(pwd)/scripts/ArrayFire-v3.9.0_Linux_x86_64.sh"
    
    if [ -f "$arrayfire_installer" ]; then
        print_status "Found ArrayFire installer in scripts/ directory"
        print_status "Installing ArrayFire to external/ArrayFire/ directory..."
        
        # Create ArrayFire directory
        mkdir -p external/ArrayFire
        cd external/ArrayFire
        
        # Run installer without --exclude-subdir to create proper subdirectory
        chmod +x "$arrayfire_installer"
        if "$arrayfire_installer"; then
            print_status "ArrayFire installer completed"
        else
            print_warning "ArrayFire installer returned non-zero exit code"
        fi
        
        cd ../..
        
        # Check for various possible directory names and installation patterns
        local possible_dirs=(
            "external/ArrayFire/ArrayFire-3.9.0-Linux"
            "external/ArrayFire-3.9.0-Linux"
            "external/ArrayFire"
            "external/arrayfire"
        )
        
        local found_dir=""
        for dir in "${possible_dirs[@]}"; do
            if [ -d "$dir" ] && ([ -f "$dir/lib64/libaf.so" ] || [ -f "$dir/lib/libaf.so" ]); then
                found_dir="$dir"
                break
            fi
        done
        
        if [ -n "$found_dir" ]; then
            export ARRAYFIRE_PATH="$found_dir"
            # Ensure lib directory exists to prevent libtool errors
            if [ ! -d "$found_dir/lib" ]; then
                mkdir -p "$found_dir/lib"
                print_status "Created missing lib directory in ArrayFire"
            fi
            print_success "ArrayFire installed successfully: $ARRAYFIRE_PATH"
        else
            print_warning "ArrayFire installation may have failed - checking what was extracted..."
            print_status "Contents of external/ directory:"
            ls -la external/ || true
        fi
    else
        print_warning "ArrayFire installer not found in scripts/ directory"
        print_status "To enable GPU acceleration:"
        print_status "  1. Download ArrayFire-v3.9.0_Linux_x86_64.sh from: https://arrayfire.com/download/"
        print_status "  2. Copy to: scripts/ArrayFire-v3.9.0_Linux_x86_64.sh"  
        print_status "  3. Re-run this script"
        print_status ""
        print_status "Note: ArrayFire installer (~2GB) is not included in the repository to keep it lightweight."
        print_status "GPU acceleration is optional - TempoNest will work without it."
    fi
}

# Function to setup PolyChord (download and build)
setup_polychord() {
    print_status "Setting up PolyChordLite for nested sampling..."
    
    local polychord_dir="$(pwd)/external/PolyChordLite"
    
    # Check if PolyChordLite already exists
    if [ -d "$polychord_dir" ] && [ -f "$polychord_dir/Makefile" ]; then
        print_status "PolyChordLite directory already exists"
        
        # Check if it's properly built
        if [ -f "$polychord_dir/lib/libchord.so" ] || [ -f "$polychord_dir/lib/libchord.a" ]; then
            print_success "PolyChordLite already built"
            export POLYCHORD_DIR="$polychord_dir"
            return 0
        else
            print_status "PolyChordLite exists but not built - will compile"
            export POLYCHORD_DIR="$polychord_dir"
            return 0
        fi
    fi
    
    # Download PolyChordLite from GitHub
    print_status "Downloading PolyChordLite from GitHub..."
    
    # Create external directory if it doesn't exist
    mkdir -p external
    
    # Clone PolyChordLite
    if command_exists git; then
        cd external
        if git clone https://github.com/PolyChord/PolyChordLite.git PolyChordLite; then
            print_success "PolyChordLite downloaded successfully"
            cd ..
            export POLYCHORD_DIR="$polychord_dir"
        else
            print_error "Failed to clone PolyChordLite repository"
            cd ..
            exit 1
        fi
    else
        print_error "Git not found - cannot download PolyChordLite"
        print_status "Please install git or manually download PolyChordLite to external/PolyChordLite/"
        exit 1
    fi
}

# Function to check required tools
check_dependencies() {
    print_status "Checking build dependencies..."
    
    local missing_deps=()
    
    if ! command_exists gcc; then
        missing_deps+=("gcc")
    fi
    
    if ! command_exists g++; then
        missing_deps+=("g++")
    fi
    
    if ! command_exists gfortran; then
        missing_deps+=("gfortran")
    fi
    
    if ! command_exists make; then
        missing_deps+=("make")
    fi
    
    if ! command_exists autoconf; then
        missing_deps+=("autoconf")
    fi
    
    if ! command_exists automake; then
        missing_deps+=("automake")
    fi
    
    if ! command_exists libtoolize; then
        missing_deps+=("libtool")
    fi
    
    # Check for BLAS library (OpenBLAS)
    if ! ldconfig -p | grep -q "libopenblas"; then
        missing_deps+=("libopenblas-dev")
    fi
    
    if [ ${#missing_deps[@]} -ne 0 ]; then
        print_error "Missing dependencies: ${missing_deps[*]}"
        print_status "Please install with: sudo apt-get install ${missing_deps[*]}"
        exit 1
    fi
    
    print_success "All build dependencies found"
}

# Function to clean previous builds
clean_build() {
    print_status "Cleaning previous build artifacts..."
    
    # Clean autotools generated files
    if [ -f "Makefile" ]; then
        make clean >/dev/null 2>&1 || true
    fi
    
    # Remove generated files (they'll be recreated)
    rm -f configure config.h config.status Makefile.in Makefile stamp-h1 libtool
    rm -f config.log config.h.in aclocal.m4
    rm -rf autom4te.cache build/config
    
    # Create build directory structure
    mkdir -p build/{obj,lib,temp}
    
    print_success "Build environment cleaned"
}

# Function to wipe everything and restore to clean source state
wipe_everything() {
    print_status "WIPING REPOSITORY - Restoring to clean source state..."
    
    # Remove all output and result files
    print_status "Removing TempoNest output files..."
    find . -name "*resume*" -delete 2>/dev/null || true
    find . -name "*TNest*" -delete 2>/dev/null || true
    find . -name "*phys_live*" -delete 2>/dev/null || true
    find . -name "*IS.points*" -delete 2>/dev/null || true
    find . -name "*ev.dat" -delete 2>/dev/null || true
    find . -name "*live.points" -delete 2>/dev/null || true
    find . -name "*post_equal_weights.dat" -delete 2>/dev/null || true
    find . -name "*stats.dat" -delete 2>/dev/null || true
    find . -name "*summary.txt" -delete 2>/dev/null || true
    
    # Remove result directories
    rm -rf results*/
    rm -rf chains/
    
    # Remove external dependencies (except bundled MultiNest)
    print_status "Removing external dependencies (preserving bundled MultiNest)..."
    if [ -d "external/ArrayFire" ]; then
        rm -rf external/ArrayFire/
    fi
    if [ -d "external/PolyChordLite" ]; then
        rm -rf external/PolyChordLite/
    fi
    
    # Remove all build artifacts
    print_status "Removing build artifacts..."
    rm -f configure config.h config.status Makefile.in Makefile stamp-h1 libtool
    rm -f config.log config.h.in aclocal.m4
    rm -rf autom4te.cache build/
    
    # Remove any compiled objects
    print_status "Removing compiled objects..."
    find . -name "*.o" -delete 2>/dev/null || true
    find . -name "*.so" -delete 2>/dev/null || true
    find . -name "*.a" -delete 2>/dev/null || true
    find . -name "*.mod" -delete 2>/dev/null || true
    find . -name "*.dirstamp" -delete 2>/dev/null || true
    
    # Remove any installed plugins
    if [ -n "$TEMPO2" ] && [ -d "$TEMPO2/plugins" ]; then
        print_status "Removing existing TempoNest plugins..."
        rm -f "$TEMPO2/plugins/temponest_"*"_plug.t2" 2>/dev/null || true
    fi
    
    # Remove generated environment file
    rm -f temponest_env.sh
    
    # Remove temporary files
    find . -name "*.tmp" -delete 2>/dev/null || true
    find . -name "*~" -delete 2>/dev/null || true
    find . -name ".DS_Store" -delete 2>/dev/null || true
    
    print_success "Repository wiped - restored to clean source state (source code + bundled MultiNest only)"
}

# Function to generate build system
generate_build_system() {
    print_status "Generating build system with autotools..."
    
    # Run autogen script
    if [ -f "scripts/autogen.sh" ]; then
        chmod +x scripts/autogen.sh
        ./scripts/autogen.sh
    else
        print_error "autogen.sh not found in scripts/ directory"
        exit 1
    fi
    
    print_success "Build system generated"
}

# Function to build bundled MultiNest
build_multinest() {
    print_status "Setting up bundled MultiNest..."
    
    if [ ! -d "external/MultiNest" ]; then
        print_error "MultiNest source not found in external/ directory - repo may be incomplete"
        exit 1
    fi
    
    # Check if bundled libraries already exist
    if [ -f "external/MultiNest/libnest3.a" ] || [ -f "external/MultiNest/libnest3.so" ]; then
        print_success "MultiNest libraries already built (libnest3)"
        return 0
    fi
    
    # Check if new format exists as fallback
    if [ -f "external/MultiNest/lib/libmultinest.a" ] || [ -f "external/MultiNest/lib/libmultinest.so" ]; then
        print_success "MultiNest libraries already built (libmultinest)"
        return 0
    fi
    
    # Need to build MultiNest from source
    print_status "Building MultiNest from source..."
    cd external/MultiNest
    
    # Check if source files exist
    if [ ! -f "nested.F90" ]; then
        print_error "MultiNest source files not found"
        cd ../..
        exit 1
    fi
    
    # Build using the existing Makefile
    if [ -f "Makefile" ]; then
        print_status "Building with existing Makefile..."
        make clean >/dev/null 2>&1 || true
        if make libnest3.so; then
            print_success "MultiNest built successfully (libnest3.so)"
            cd ../..
            return 0
        elif make libnest3.a; then
            print_success "MultiNest built successfully (libnest3.a)"
            cd ../..
            return 0
        else
            print_error "MultiNest build failed"
            cd ../..
            exit 1
        fi
    else
        print_error "No Makefile found in MultiNest directory"
        cd ../..
        exit 1
    fi
}

# Function to build PolyChordLite
build_polychord() {
    print_status "Building PolyChordLite..."
    
    if [ ! -d "external/PolyChordLite" ]; then
        print_error "PolyChordLite source not found in external/ directory"
        exit 1
    fi
    
    # Check if already built
    if [ -f "external/PolyChordLite/lib/libchord.so" ] || [ -f "external/PolyChordLite/lib/libchord.a" ]; then
        print_success "PolyChordLite already built"
        return 0
    fi
    
    # Build PolyChordLite
    print_status "Building PolyChordLite from source..."
    cd external/PolyChordLite
    
    # Check if Makefile exists
    if [ ! -f "Makefile" ]; then
        print_error "PolyChordLite Makefile not found"
        cd ../..
        exit 1
    fi
    
    # Clean previous builds
    make clean >/dev/null 2>&1 || true
    
    # Build the library
    if make; then
        print_success "PolyChordLite built successfully"
        cd ../..
        return 0
    else
        print_error "PolyChordLite build failed"
        cd ../..
        exit 1
    fi
}

# Function to setup library paths
setup_library_paths() {
    print_status "Setting up library paths..."
    
    # Ensure MULTINEST_DIR is set
    if [ -z "$MULTINEST_DIR" ]; then
        export MULTINEST_DIR="$(pwd)/external/MultiNest"
        print_status "Set MULTINEST_DIR to: $MULTINEST_DIR"
    fi
    
    # Ensure POLYCHORD_DIR is set if PolyChord exists
    if [ -z "$POLYCHORD_DIR" ]; then
        local polychord_dir="$(pwd)/external/PolyChordLite"
        if [ -d "$polychord_dir" ]; then
            export POLYCHORD_DIR="$polychord_dir"
            print_status "Set POLYCHORD_DIR to: $POLYCHORD_DIR"
        fi
    fi
    
    # Ensure ARRAYFIRE_PATH is set if ArrayFire exists
    if [ -z "$ARRAYFIRE_PATH" ]; then
        local possible_arrayfire_dirs=(
            "$(pwd)/external/ArrayFire"
            "$(pwd)/external/ArrayFire-3.9.0-Linux"
        )
        
        for dir in "${possible_arrayfire_dirs[@]}"; do
            if [ -d "$dir" ] && ([ -f "$dir/lib64/libaf.so" ] || [ -f "$dir/lib/libaf.so" ]); then
                export ARRAYFIRE_PATH="$dir"
                print_status "Set ARRAYFIRE_PATH to: $ARRAYFIRE_PATH"
                break
            fi
        done
    fi
    
    # Set up LDFLAGS and LD_LIBRARY_PATH for MultiNest
    local multinest_path="$MULTINEST_DIR"
    
    # Check if bundled MultiNest (libraries in root directory)
    if [ -f "$multinest_path/libnest3.so" ] || [ -f "$multinest_path/libnest3.a" ]; then
        # Bundled MultiNest structure with libraries in root
        export LDFLAGS="-L$multinest_path $LDFLAGS"
        export LD_LIBRARY_PATH="$multinest_path:$LD_LIBRARY_PATH"
    elif [ -d "$multinest_path/lib" ]; then
        # New MultiNest structure with lib/ subdirectory
        export LDFLAGS="-L$multinest_path/lib $LDFLAGS"
        export LD_LIBRARY_PATH="$multinest_path/lib:$LD_LIBRARY_PATH"
    else
        # Fallback to root directory
        export LDFLAGS="-L$multinest_path $LDFLAGS"
        export LD_LIBRARY_PATH="$multinest_path:$LD_LIBRARY_PATH"
    fi
    
    # Add PolyChord paths if available
    if [ -n "$POLYCHORD_DIR" ]; then
        export LDFLAGS="-L$POLYCHORD_DIR/lib $LDFLAGS"
        export LD_LIBRARY_PATH="$POLYCHORD_DIR/lib:$LD_LIBRARY_PATH"
    fi
    
    # Add ArrayFire paths if available
    if [ -n "$ARRAYFIRE_PATH" ]; then
        export LDFLAGS="-L$ARRAYFIRE_PATH/lib64 -L$ARRAYFIRE_PATH/lib $LDFLAGS"
        export LD_LIBRARY_PATH="$ARRAYFIRE_PATH/lib64:$ARRAYFIRE_PATH/lib:$LD_LIBRARY_PATH"
    fi
    
    print_success "Library paths configured"
    print_status "MULTINEST_DIR: $MULTINEST_DIR"
    if [ -n "$POLYCHORD_DIR" ]; then
        print_status "POLYCHORD_DIR: $POLYCHORD_DIR"
    fi
    if [ -n "$ARRAYFIRE_PATH" ]; then
        print_status "ARRAYFIRE_PATH: $ARRAYFIRE_PATH"
    fi
}

# Function to configure build
configure_build() {
    print_status "Configuring build..."
    
    local configure_args=""
    
    # Add BLAS configuration (from your notes)
    configure_args="$configure_args --with-blas=openblas"
    
    # Add ArrayFire if available
    if [ -n "$ARRAYFIRE_PATH" ]; then
        configure_args="$configure_args --with-arrayfire=$ARRAYFIRE_PATH"
        print_status "Configuring with ArrayFire support"
    fi
    
    # Add MultiNest
    if [ -n "$MULTINEST_DIR" ]; then
        configure_args="$configure_args --with-multinest=$MULTINEST_DIR"
    fi
    
    # Add PolyChord if available
    if [ -n "$POLYCHORD_DIR" ]; then
        configure_args="$configure_args --with-polychord=$POLYCHORD_DIR"
        print_status "Configuring with PolyChord support"
    fi
    
    # Ensure environment variables are exported for configure
    if [ -n "$MULTINEST_DIR" ]; then
        export MULTINEST_DIR
        print_status "Exported MULTINEST_DIR=$MULTINEST_DIR"
    fi
    
    if [ -n "$POLYCHORD_DIR" ]; then
        export POLYCHORD_DIR
        print_status "Exported POLYCHORD_DIR=$POLYCHORD_DIR"
    fi
    
    if [ -n "$ARRAYFIRE_PATH" ]; then
        export ARRAYFIRE_PATH
        print_status "Exported ARRAYFIRE_PATH=$ARRAYFIRE_PATH"
    fi
    
    # Run configure
    print_status "Running: ./configure $configure_args"
    ./configure $configure_args
    
    print_success "Configuration complete"
}

# Function to compile TempoNest
compile_temponest() {
    print_status "Compiling TempoNest..."
    
    # Get number of CPU cores for parallel compilation
    local num_cores=$(nproc 2>/dev/null || echo 4)
    
    print_status "Using $num_cores parallel jobs"
    make -j$num_cores temponest
    
    print_success "TempoNest compilation complete"
}

# Function to run tests
run_tests() {
    print_status "Running tests..."
    
    if make test >/dev/null 2>&1; then
        print_success "All tests passed"
    else
        print_warning "Some tests failed - check test output"
    fi
}

# Function to install plugin
install_plugin() {
    print_status "Installing TempoNest plugin to Tempo2..."
    
    make temponest-install
    
    # Verify installation
    local plugin_dir="$TEMPO2/plugins"
    if [ -f "$plugin_dir/temponest_"*"_plug.t2" ]; then
        print_success "Plugin installed successfully to $plugin_dir"
    else
        print_error "Plugin installation failed"
        exit 1
    fi
}

# Function to verify installation
verify_installation() {
    print_status "Verifying installation..."
    
    # Test if tempo2 can load the plugin
    local test_output
    test_output=$(echo "quit" | timeout 10s $TEMPO2/bin/tempo2 -gr temponest 2>&1 || true)
    
    if echo "$test_output" | grep -q "TempoNest"; then
        print_success "TempoNest plugin loaded successfully"
    else
        print_warning "Could not verify plugin loading - but installation may still be successful"
    fi
}

# Function to create environment script for runtime
create_environment_script() {
    print_status "Creating comprehensive environment script..."
    
    local env_file="temponest_env.sh"
    
    cat > "$env_file" << EOF
#!/bin/bash
# TempoNest Environment Setup
# Generated by build_temponest.sh on $(date)
# Source this file before using TempoNest: source ./temponest_env.sh

# ===== TempoNest Dependencies =====
EOF
    
    # MultiNest - use absolute paths (bundled with repo)
    if [ -n "$MULTINEST_DIR" ] && [ -d "$MULTINEST_DIR" ]; then
        local abs_multinest_dir=$(cd "$MULTINEST_DIR" && pwd)
        echo "export MULTINEST_DIR=\"$abs_multinest_dir\"" >> "$env_file"
        
        # Check if bundled MultiNest (libraries in root directory)
        if [ -f "$abs_multinest_dir/libnest3.so" ] || [ -f "$abs_multinest_dir/libnest3.a" ]; then
            echo "export LD_LIBRARY_PATH=\"$abs_multinest_dir:\$LD_LIBRARY_PATH\"" >> "$env_file"
        elif [ -d "$abs_multinest_dir/lib" ]; then
            echo "export LD_LIBRARY_PATH=\"$abs_multinest_dir/lib:\$LD_LIBRARY_PATH\"" >> "$env_file"
        else
            echo "export LD_LIBRARY_PATH=\"$abs_multinest_dir:\$LD_LIBRARY_PATH\"" >> "$env_file"
        fi
    fi
    
    # PolyChord - use absolute paths
    if [ -n "$POLYCHORD_DIR" ] && [ -d "$POLYCHORD_DIR" ]; then
        local abs_polychord_dir=$(cd "$POLYCHORD_DIR" && pwd)
        echo "export POLYCHORD_DIR=\"$abs_polychord_dir\"" >> "$env_file"
        echo "export LD_LIBRARY_PATH=\"$abs_polychord_dir/lib:\$LD_LIBRARY_PATH\"" >> "$env_file"
    fi
    
    # ArrayFire - use absolute paths
    if [ -n "$ARRAYFIRE_PATH" ] && [ -d "$ARRAYFIRE_PATH" ]; then
        local abs_arrayfire_path=$(cd "$ARRAYFIRE_PATH" && pwd)
        echo "export ARRAYFIRE_PATH=\"$abs_arrayfire_path\"" >> "$env_file"
        echo "export LD_LIBRARY_PATH=\"$abs_arrayfire_path/lib64:$abs_arrayfire_path/lib:\$LD_LIBRARY_PATH\"" >> "$env_file"
    fi
    
    # Tempo2 (detect if available)
    if [ -n "$TEMPO2" ] && [ -d "$TEMPO2" ]; then
        echo "export TEMPO2=\"$TEMPO2\"" >> "$env_file"
        echo "export PATH=\"$TEMPO2/bin:\$PATH\"" >> "$env_file"
        echo "export LD_LIBRARY_PATH=\"$TEMPO2/lib:\$LD_LIBRARY_PATH\"" >> "$env_file"
    elif command -v tempo2 >/dev/null 2>&1; then
        local tempo2_bin=$(which tempo2)
        local tempo2_root=$(dirname $(dirname $tempo2_bin))
        echo "export TEMPO2=\"$tempo2_root\"" >> "$env_file"
        echo "export PATH=\"$tempo2_root/bin:\$PATH\"" >> "$env_file"
        if [ -d "$tempo2_root/lib" ]; then
            echo "export LD_LIBRARY_PATH=\"$tempo2_root/lib:\$LD_LIBRARY_PATH\"" >> "$env_file"
        fi
    fi
    
    # Add any additional library paths that might be needed
    cat >> "$env_file" << 'EOF'

# Remove duplicate entries from LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | awk '!seen[$0]++' | tr '\n' ':' | sed 's/:$//')

# ===== Status Messages =====
echo "==============================================="
echo "TempoNest Environment Loaded"
echo "==============================================="
echo "TEMPO2: $TEMPO2"
echo "MULTINEST_DIR: $MULTINEST_DIR"
echo "POLYCHORD_DIR: $POLYCHORD_DIR"
echo "ARRAYFIRE_PATH: $ARRAYFIRE_PATH"
echo "LD_LIBRARY_PATH configured for runtime"
echo "==============================================="
echo "You can now run: tempo2 -gr temponest ..."
echo "==============================================="
EOF
    
    chmod +x "$env_file"
    print_success "Comprehensive environment script created: $env_file"
    print_status "Users should run: source ./$env_file"
}

# Function to display usage information
usage() {
    cat << EOF
TempoNest Build Script

Usage: $0 [options]

Options:
    --clean         Clean build directory before building (keeps compiled artifacts)
    --wipe          Restore repository to clean source state (remove all build artifacts)
    --no-tests      Skip running tests
    --no-install    Skip plugin installation
    --help          Show this help message

Dependencies:
    - Uses bundled MultiNest (included in repository)
    - Downloads PolyChordLite from GitHub if not present
    - Uses ArrayFire installer from scripts/ directory if available
    - Auto-detects Tempo2 installation

Example:
    $0 --wipe       # Complete fresh build (clean source state)
    $0 --clean      # Incremental build (keeps compiled artifacts)
    $0              # Build with existing configuration

EOF
}

# Parse command line arguments
CLEAN_BUILD=false
WIPE_EVERYTHING=false
RUN_TESTS=true
INSTALL_PLUGIN=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --wipe)
            WIPE_EVERYTHING=true
            shift
            ;;
        --no-tests)
            RUN_TESTS=false
            shift
            ;;
        --no-install)
            INSTALL_PLUGIN=false
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Main build process
main() {
    print_status "Starting TempoNest build process..."
    print_status "Build directory: $(pwd)"
    
    # Wipe everything if requested and exit
    if [ "$WIPE_EVERYTHING" = true ]; then
        wipe_everything
        print_status ""
        print_status "Repository wiped successfully!"
        print_status "To build TempoNest, run: ./scripts/build_temponest.sh"
        return 0
    fi
    
    # Initialize clean build environment (autonomous operation)
    initialize_build_environment
    
    # Check prerequisites
    check_dependencies
    setup_dependencies
    
    # Clean if requested
    if [ "$CLEAN_BUILD" = true ]; then
        clean_build
    fi
    
    # Build process
    build_multinest
    build_polychord
    setup_library_paths
    generate_build_system
    configure_build
    compile_temponest
    
    # Optional steps
    if [ "$RUN_TESTS" = true ]; then
        run_tests
    fi
    
    if [ "$INSTALL_PLUGIN" = true ]; then
        install_plugin
        verify_installation
    fi
    
    print_success "TempoNest build complete!"
    
    # Create environment script for future use
    create_environment_script
    
    print_status ""
    print_status "You can now use TempoNest with:"
    print_status "  tempo2 -gr temponest -f par.file tim.file -Cfile config.json"
    print_status ""
    print_status "If you encounter library loading issues, run:"
    print_status "  source ./temponest_env.sh"
    print_status ""
    
    if [ -n "$ARRAYFIRE_PATH" ]; then
        print_status "GPU acceleration is enabled via ArrayFire"
    else
        print_status "GPU acceleration is disabled (ArrayFire not found)"
    fi
}

# Run main function
main "$@"