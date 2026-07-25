@echo off
echo ==============================
echo Building Parallel Matrix Inverse
echo ==============================

g++ -fopenmp ^
main.cpp ^
Metrices\*.cpp ^
-o main.exe


if %ERRORLEVEL% neq 0 (
    echo.
    echo Build FAILED!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo Build SUCCESS!
echo.

main.exe

echo.
pause