@echo off
rem ============================================================
rem  build.bat              : debug build (bounds check, backtrace)
rem  build.bat release      : optimized build
rem  build.bat release omp  : optimized + OpenMP parallel
rem  build.bat omp          : debug + OpenMP
rem
rem  OpenMP threads: $env:OMP_NUM_THREADS = 8  (PowerShell, before running)
rem
rem  SRC order = dependency order
rem  (a module must come before any file that uses it)
rem ============================================================

set SRC=src/euler.f90 src/case_setup.f90 src/grid.f90 src/flow.f90 src/solver.f90 src/monitor.f90 src/output.f90 src/main.f90

set FLAGS=-O0 -g -Wall -fcheck=all -fbacktrace
set OMP=
for %%a in (%*) do (
    if "%%a"=="release" set FLAGS=-O3
    if "%%a"=="omp" set OMP=-fopenmp
)

if not exist build mkdir build
gfortran %FLAGS% %OMP% -Jbuild %SRC% -o main.exe
if errorlevel 1 (
    echo === BUILD FAILED ===
    exit /b 1
)
echo === BUILD OK: main.exe  [%FLAGS% %OMP%] ===
