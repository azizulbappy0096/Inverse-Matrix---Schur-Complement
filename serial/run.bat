@echo off
echo ==============================
echo Building Serial Matrix Inverse
echo ==============================

g++ main.cpp ^
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