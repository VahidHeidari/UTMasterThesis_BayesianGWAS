@echo off

echo Clean from root directory . . .
del /s *~ *.pyc *.gtp

echo Call ReverseMapTest clean script . . .
cd ReverseMapTest
call clean.bat
cd ..

echo Done!
