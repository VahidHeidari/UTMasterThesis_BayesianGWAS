#!/bin/bash

ARCHIVE_FILE=DataGen-`date +%Y%m%d-%a`.tar.gz
tar -czvf $ARCHIVE_FILE *.h *.cpp *.bat *.sh *.sln *.vcxproj *.filters simulation_parameters.txt
sleep 3

