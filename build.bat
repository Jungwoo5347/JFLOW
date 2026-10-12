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

set SRC=euler.f90 case_setup.f90 grid.f90 flow.f90 solver.f90 output.f90 main.f90

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
