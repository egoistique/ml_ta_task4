@echo off
set CD=%~dp0

::set javaDir="C:\Program Files\Java\jdk-16.0.2\bin"

set PATH=%javaDir%;%PATH%

set cp=%cd%\out\production\Arr2JTableDemo

set main=ru.vsu.cs.course1.Program 


java -classpath %cp% %main% %*
