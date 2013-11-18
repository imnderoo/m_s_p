@echo off
echo.
REM Ping provides a delay
ping -n 2 127.0.0.1 > nul:

C:\cygwin\bin\bash.exe --login -c "/usr/local/bin/dos2cygwin.sh"

ping -n 2 127.0.0.1 > nul:
echo.
echo. Create Report Script completed.
echo. Press any key to exit...
pause > null
