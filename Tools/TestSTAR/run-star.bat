@echo off

set FIRST_DIR=%PATH:~0,15%
set CYGWIN_PATH=C:\cygwin64\bin
set STAR_EXE=E:\C++\STAR\source\STAR.exe

if %FIRST_DIR% NEQ %CYGWIN_PATH% ( call :SET_NEW_PATH )

%STAR_EXE% %*
goto:EOF



:SET_NEW_PATH
	echo Cygwin bin directory is added to %%PATH%% environment variable.
	set PATH=%CYGWIN_PATH%;%PATH%
goto:EOF

