@echo off
rem ============================================================
rem  build.bat              : debug build (bounds check, backtrace)
rem  build.bat release      : optimized build
rem  build.bat release omp  : optimized + OpenMP parallel
rem  build.bat omp          : debug + OpenMP
rem  build.bat 3d ...       : 3D solver (src3d/) -> main3d.exe, same options
rem                           e.g.  build.bat 3d release omp
rem
rem  OpenMP threads: $env:OMP_NUM_THREADS = 8  (PowerShell, before running)
rem
rem  SRC order = dependency order
rem  (a module must come before any file that uses it)
rem ============================================================

set SRC=src/euler.f90 src/case_setup.f90 src/grid.f90 src/flow.f90 src/solver.f90 src/monitor.f90 src/output.f90 src/main.f90
set SRC3D=src3d/euler.f90 src3d/case_setup.f90 src3d/grid.f90 src3d/flow.f90 src3d/solver.f90 src3d/implicit.f90 src3d/monitor.f90 src3d/output.f90 src3d/main.f90

set FLAGS=-O0 -g -Wall -fcheck=all -fbacktrace
set OMP=
set EXE=main.exe
set MODDIR=build
for %%a in (%*) do (
    if "%%a"=="release" set FLAGS=-O3
    if "%%a"=="omp" set OMP=-fopenmp
    if "%%a"=="3d" (
        set SRC=%SRC3D%
        set EXE=main3d.exe
        set MODDIR=build3d
    )
)

if not exist %MODDIR% mkdir %MODDIR%
gfortran %FLAGS% %OMP% -J%MODDIR% %SRC% -o %EXE%
if errorlevel 1 (
    echo === BUILD FAILED ===
    exit /b 1
)
echo === BUILD OK: %EXE%  [%FLAGS% %OMP%] ===
