@echo off

set Rscript=%1%
set R_HOME=%2%
set WDrive=%3%
set WDir=%4%
set Master=%5%
set MapDrive=%6%
set RemotePath=%7%

if %MapDrive%==FALSE goto last 
net use %WDrive%: %RemotePath%
if not errorlevel 0 goto last

cd /D %WDir%
%R_HOME%\bin\x64\Rscript.exe --vanilla %Rscript%
goto :EOF

:last
if %Master% == %COMPUTERNAME% (
cd /D %WDir%
if not errorlevel 0 cd /D e:\
%R_HOME%\bin\x64\Rscript.exe --vanilla %Rscript%
goto:eof 
)

:home
cd /D e:\
%R_HOME%\bin\x64\Rscript.exe --vanilla %Rscript%
