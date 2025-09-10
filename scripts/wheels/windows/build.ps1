param([System.String]$cotire="OFF")

Write-Host "building kratos version: $env:KRATOS_VERSION"

$pythonVersion = $env:PYTHON_VERSION
$python = $pythonVersion.Replace(".", "")

$kratosRoot = Resolve-Path "."
$env:kratos_root = $kratosRoot
$wheelRoot = Join-Path $kratosRoot "wheel"
$wheelOutDir = Join-Path $env:GITHUB_WORKSPACE "data_swap_guest"
mkdir "$($wheelOutDir)"

$numcores = (Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors

function exec_build($python, $pythonPath) {
    cmd.exe /c "call configure_wheel.bat $($pythonPath) $($kratosRoot) OFF"
    cmake --build "$($kratosRoot)/build/Release" --target install -- /property:configuration=Release /p:Platform=x64
}

function exec_build_cotire($python, $pythonPath) {
    cmd.exe /c "call configure_wheel.bat $($pythonPath) $($kratosRoot) ON"
    cmake --build "$($kratosRoot)/build/Release" --target install -- /property:configuration=Release /p:Platform=x64
}

function  setup_wheel_dir {
    cd $kratosRoot
    mkdir "$($wheelRoot)"
    cp scripts\wheels\setup.py "$($wheelRoot)\setup.py"
    mkdir "$($wheelRoot)\KratosMultiphysics"
    mkdir "$($wheelRoot)\KratosMultiphysics\.libs"
}

function build_core_wheel ($pythonLocation, $prefixLocation) {
    setup_wheel_dir
    cd $kratosRoot

    cp $prefixLocation\KratosMultiphysics\*         "$($wheelRoot)\KratosMultiphysics"
    cp $kratosRoot\kratos\KratosMultiphysics.json   "$($wheelRoot)\wheel.json"

    cd $wheelRoot

    & $pythonLocation setup.py bdist_wheel

    cp "$($wheelRoot)\dist\*" $wheelOutDir
    
    cd c:\
    rm -r $wheelRoot
}

function build_application_wheel ($pythonPath, $app) {
    setup_wheel_dir
    cp "$($kratosRoot)\applications\$($app)\$($app).json" "$($wheelRoot)\wheel.json"
    cd $wheelRoot
    & $pythonPath setup.py bdist_wheel #pythonpath
    cp "$($wheelRoot)\dist\*" $wheelOutDir
    cd c:\
    rm -r $wheelRoot
}

function build_core ($pythonLocation, $prefixLocation) {
    cd $kratosRoot

    Write-Host $kratosRoot
    Write-Host "$($kratosRoot)\scripts\wheels\windows\configure.bat"

    cp "$($kratosRoot)\scripts\wheels\windows\configure.bat" .\configure_wheel.bat

    Write-Host "Debugging: begin cmd.exe call"
    
    cmd.exe /c "call configure.bat $($pythonLocation) $($kratosRoot) $($prefixLocation) $($numcores)"
    cmake --build "$($kratosRoot)/build/Release" --target KratosKernel -- /property:configuration=Release /p:Platform=x64 /p:CL_MPCount=24 /m:1
}

function build_interface ($pythonLocation, $pythonPath) {
    cd $kratosRoot

    cp "$($kratosRoot)\scripts\wheels\windows\configure.bat" .\configure_wheel.bat

    cmd.exe /c "call configure_wheel.bat $($pythonLocation) $($kratosRoot) $($prefixLocation) $($numcores)"
    cmake --build "$($kratosRoot)/build/Release" --target KratosPythonInterface -- /property:configuration=Release /p:Platform=x64 /p:CL_MPCount=24 /m:1
    cmake --build "$($kratosRoot)/build/Release" --target install -- /property:configuration=Release /p:Platform=x64 /p:CL_MPCount=24 /m:1
}

# Core can be built independently of the python version.
# Install path should be useless here.
# Write-Host "Starting core build"
# # build_core "$($env:pythonRoot)\39\python.exe" ${KRATOS_ROOT}/bin/core
# Write-Host "Finished core build"


Write-Host "Beginning build for python $($python)"
$env:python = $python

cd $kratosRoot

$pythonLocation = (Get-Command python).Source
$prefixLocation = "$($kratosRoot)\bin\Release\python_$($python)"


build_core $pythonLocation $prefixLocation
build_interface $pythonLocation $prefixLocation

Write-Host "Finished build"

Write-Host "Building Core Wheel"
build_core_wheel $pythonLocation $prefixLocation

$applications = Get-ChildItem -Path "${prefixLocation}\applications" -Directory -Force -ErrorAction SilentlyContinue

Write-Host "Building App Wheels"
foreach($app in $applications) {
    build_application_wheel $pythonLocation $app.Name
}

Write-Host "Finished wheel construction for python $($python)"

