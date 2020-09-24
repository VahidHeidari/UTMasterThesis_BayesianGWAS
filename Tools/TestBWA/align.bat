@echo off

python mk_samples.py					^
	&& mk-index.bat						^
	&& mk-map.bat						^
	&& python align_cover.py

