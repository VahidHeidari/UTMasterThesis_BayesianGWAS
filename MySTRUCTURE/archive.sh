#!/bin/bash

ARCHIVE_NAME=MySTRUCTURE-`date +%Y%m%d-%a`.tar.gz
tar --exclude=Debug --exclude=Release -cvzf $ARCHIVE_NAME		\
	*.bat *.sh .gitignore *.sln *.vcxproj *.filters	config-help.txt	\
	Libs NoAdmixture Admixture ExhMotahari FastSTRUCTURE

echo Done . . .
sleep 3

