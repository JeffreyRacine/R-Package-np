@echo off

set Rscript=%1%
set R_HOME=%2%
set WDrive=%3%
set WDir=%4%
set Master=%5%
set MapDrive=%6%
set RemotePath=%7%

if %MapDrive%==FALSE goto last 
set FirstOne=FALSE
if exist %WDrive%:\NUL goto next
set FirstOne=TRUE
net use %WDrive%: %RemotePath%
if not errorlevel 0 set FirstOne=FALSE

:next
cd /D %WDir%
%R_HOME%\bin\i386\Rscript.exe --vanilla %Rscript%

if %FirstOne%==FALSE goto end 
net use /delete %WDrive%:
goto:eof

:last
if %Master% == %COMPUTERNAME% (
%R_HOME%\bin\i386\Rscript.exe --vanilla %Rscript%
goto:eof 
)

cd /D %TMP%
%R_HOME%\bin\i386\Rscript.exe --vanilla %Rscript%

:end