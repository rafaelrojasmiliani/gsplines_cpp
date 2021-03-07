
@echo off
for /F "tokens=1,2,3" %%A in ('netsh interface ip show addresses "vEthernet (Default Switch)" ^| findstr "IP Address"') DO (set xserverip=%%C)

rem START /B C:\"Program Files"\VcXsrv\vcxsrv.exe
docker run -it ^
        --volume %CD%\..:/gsplinespp ^
        --user %USERNAME%:%USERNAME% ^
        --env="DISPLAY=%xserverip%:0.0" ^
        "gsplinespp" bash
