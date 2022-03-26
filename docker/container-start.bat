@echo off
for /F "tokens=1,2,3" %%A in ('netsh interface ip show addresses "vEthernet (Default Switch)" ^| findstr "IP Address"') DO (set xserverip=%%C)

rem START /B C:\"Program Files"\VcXsrv\vcxsrv.exe
docker pull rafa606/cpp-vim
docker run -it ^
        --volume %CD%\..:/workspace ^
        --env="DISPLAY=%xserverip%:0.0" ^
        --entrypoint="/bin/bash" ^
        "rafa606/cpp-vim" -c "addgroup --gid 11021 %USERNAME% --force-badname;  adduser --gecos \"\" --disabled-password  --uid 11021 --gid 11021 %USERNAME% --force-badname --home /home/%USERNAME%; install -d -m 0755 -o %USERNAME% -g 11021 /home/%USERNAME%;  usermod -a -G video %USERNAME%; echo %USERNAME% ALL=\(ALL\) NOPASSWD:ALL >> /etc/sudoers; sudo -EHu %USERNAME%  bash"
