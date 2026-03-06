#!/bin/bash
export PYTHON=$PYTHON_VERSION

echo "All envs:"
env | sort

BASE_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
export KRATOS_ROOT="$(pwd)"
WHEEL_ROOT="$KRATOS_ROOT/wheel"
WHEEL_OUT="$GITHUB_WORKSPACE/data_swap_guest"
mkdir -p $WHEEL_OUT
CORE_LIB_DIR="$KRATOS_ROOT/coreLibs"

# Created the wheel building directory.
setup_wheel_dir () {
    cd $KRATOS_ROOT
    mkdir $WHEEL_ROOT
    cp scripts/wheels/setup.py ${WHEEL_ROOT}/setup.py
    mkdir ${WHEEL_ROOT}/KratosMultiphysics
    mkdir ${WHEEL_ROOT}/KratosMultiphysics/.libs        # This dir is necessary mainly to store the shared libs in windows
}

# Creates the wheel for KratosCore.
build_core_wheel () {
    setup_wheel_dir
    cd $KRATOS_ROOT

    PREFIX_LOCATION=$1

#    mkdir ${WHEEL_ROOT}/KratosMultiphysics

    cp ${PREFIX_LOCATION}/KratosMultiphysics/*       ${WHEEL_ROOT}/KratosMultiphysics
    cp ${KRATOS_ROOT}/kratos/KratosMultiphysics.json ${WHEEL_ROOT}/wheel.json

    cd $WHEEL_ROOT

    $PYTHON_LOCATION setup.py bdist_wheel

    cd ${WHEEL_ROOT}/dist

    auditwheel repair *.whl

    mkdir $CORE_LIB_DIR
    unzip -j wheelhouse/KratosMultiphysics* 'KratosMultiphysics.libs/*' -d $CORE_LIB_DIR

    cp wheelhouse/* ${WHEEL_OUT}/
    cd
    rm -r $WHEEL_ROOT
}

# Creates the wheel for each KratosApplication.
build_application_wheel () {
    setup_wheel_dir

    cp ${KRATOS_ROOT}/applications/${1}/${1}.json ${WHEEL_ROOT}/wheel.json
    cd $WHEEL_ROOT

    $PYTHON_LOCATION setup.py bdist_wheel

    auditwheel repair dist/*.whl

    optimize_wheel

    cp ${WHEEL_ROOT}/wheelhouse/* ${WHEEL_OUT}/

    cd
    rm -r $WHEEL_ROOT
}

# Kreates the wheel bundle.
build_kratos_all_wheel () {
    setup_wheel_dir
    cp ${KRATOS_ROOT}/kratos/KratosMultiphysics-all.json ${WHEEL_ROOT}/wheel.json
    cp ${KRATOS_ROOT}/scripts/wheels/linux/setup_kratos_all.py ${WHEEL_ROOT}/setup.py
    cd ${WHEEL_ROOT}
    $PYTHON_LOCATION setup.py bdist_wheel
    cp dist/* ${WHEEL_OUT}/

    cd
    rm -r $WHEEL_ROOT
}

# Removes duplicated libraries from existing wheels.
optimize_wheel() {
    local wheel_dir="${WHEEL_ROOT}/wheelhouse"
    local tmp_dir="tmp_wheel_extract"
    local archive_name

    cd "${wheel_dir}" || return 1
    archive_name=$(ls *.whl 2>/dev/null | head -n 1)
    
    if [ -z "${archive_name}" ]; then
        echo "Error: No wheel file found in ${wheel_dir}"
        return 1
    fi

    echo "Optimizing wheel: ${archive_name}"

    mkdir -p "${tmp_dir}"
    unzip -q "${archive_name}" -d "${tmp_dir}"
    rm -f "${archive_name}"

    # 1. ROBUST FOLDER SEARCH: Find the .libs folder, whatever auditwheel named it
    local libs_dir
    libs_dir=$(find "${tmp_dir}" -type d -name "*.libs" | head -n 1)

    if [ -n "${libs_dir}" ] && [ -d "${libs_dir}" ]; then
        echo "Found auditwheel libs directory: ${libs_dir}"
        
        for library_path in "${libs_dir}"/*; do
            [ -e "${library_path}" ] || continue 
            
            local library=$(basename "${library_path}")
            
            # 2. STRIP THE AUDITWHEEL HASH: Extract 'libKratosCore' from 'libKratosCore-1a2b3c4d.so'
            local lib_prefix=$(echo "${library}" | awk -F'-' '{print $1}')

            # 3. WILDCARD CHECK: Does ANY file starting with this prefix exist in CORE_LIB_DIR?
            # We also check the excluded.txt just in case.
            if ls "${CORE_LIB_DIR}/${lib_prefix}"* 1> /dev/null 2>&1 || \
               grep -q "^${lib_prefix}" "${WHEEL_ROOT}/excluded.txt" 2>/dev/null; then
                
                echo "-- Removing ${library} (Matched core prefix: ${lib_prefix})"
                rm -f "${library_path}"

                # Safely remove from the RECORD file
                local safe_lib_name=$(echo "${library}" | sed 's/\./\\./g')
                sed -i "/${safe_lib_name}/d" "${tmp_dir}"/*.dist-info/RECORD
            else
                echo "-- Keeping ${library}"
            fi
        done
    else
        echo "-- No .libs directory found. Auditwheel may not have bundled dependencies here."
    fi

    # Clean the possible copies done in setup.py
    rm -f "${tmp_dir}/KratosMultiphysics/.libs/libKratos"*

    # Re-zip the wheel safely
    cd "${tmp_dir}" || return 1
    zip -qr "../${archive_name}" .
    cd ..
    
    rm -rf "${tmp_dir}"
    echo "Optimization complete: ${archive_name}"
}

# Buils the KratosXCore components for the kernel and applications
build_core () {
    PYTHON_LOCATION=$1
    PREFIX_LOCATION=$2

	cp $KRATOS_ROOT/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION $PREFIX_LOCATION

    cmake --build "${KRATOS_ROOT}/build/Release" --target KratosKernel -- -j$(nproc)
}

# Buils the KratosXInterface components for the kernel and applications given an specific version of python
build_interface () {

    PYTHON_LOCATION=$1
    PREFIX_LOCATION=$2

	cp $KRATOS_ROOT/scripts/wheels/linux/configure.sh ./configure.sh
	chmod +x configure.sh
	./configure.sh $PYTHON_LOCATION $PREFIX_LOCATION

    cmake --build "${KRATOS_ROOT}/build/Release" --target KratosPythonInterface -- -j$(nproc)
    cmake --build "${KRATOS_ROOT}/build/Release" --target install -- -j$(nproc)
}

# Core can be build independently of the python version.
# Install path should be useless here.


PYTHON_LOCATION=$(which python$PYTHON_VERSION)

$PYTHON_LOCATION -m pip install setuptools wheel auditwheel

echo "Starting core build"
build_core $PYTHON_LOCATION ${KRATOS_ROOT}/bin/core
echo "Finished core build"


  echo "Starting build for python${PYTHON_VERSION}"
  PREFIX_LOCATION="$KRATOS_ROOT/bin/Release/python_$PYTHON_VERSION"

  $PYTHON_LOCATION -m pip install mypy

  build_interface $PYTHON_LOCATION $PREFIX_LOCATION

cd $KRATOS_ROOT
export HASH=$(git show -s --format=%h) # Used in version number
export LD_LIBRARY_PATH=${PREFIX_LOCATION}/libs:$BASE_LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

  echo "Building Core Wheel"
  build_core_wheel $PREFIX_LOCATION

  echo "Building App Wheels"
  for APPLICATION in $(ls -d ${PREFIX_LOCATION}/applications/*)
  do
      APPNAME=$(basename "$APPLICATION")
      echo "Building ${APPNAME} Wheel"
      build_application_wheel $APPNAME
  done

  echo "Building Bundle Wheel"
  build_kratos_all_wheel $PREFIX_LOCATION

echo finished build for python${PYTHON_VERSION}

export LD_LIBRARY_PATH=$BASE_LD_LIBRARY_PATH
