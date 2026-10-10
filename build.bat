@echo off
rem ============================================================
rem  build.bat          : debug build (bounds check, backtrace)
rem  build.bat release  : optimized build
rem
rem  SRC order = dependency order
rem  (a module must come before any file that uses it)
rem ============================================================

set SRC=param.f90 euler.f90 grid.f90 init_cond.f90 flow.f90 main.f90

set FLAGS=-O0 -g -Wall -fcheck=all -fbacktrace
if "%1"=="release" set FLAGS=-O3

if not exist build mkdir build
gfortran %FLAGS% -Jbuild %SRC% -o main.exe
if errorlevel 1 (
    echo === BUILD FAILED ===
    exit /b 1
)
echo === BUILD OK: main.exe ===
