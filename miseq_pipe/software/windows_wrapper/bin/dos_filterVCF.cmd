@echo off
echo.

C:\cygwin\bin\bash.exe --login -c "/usr/local/bin/dos_filterVCF.sh"

REM Ping provides a delay
ping -n 2 127.0.0.1 > nul:
echo.Job completed. Exiting...
ping -n 20 127.0.0.1 > nul:
echo. 