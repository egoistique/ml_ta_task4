@echo off

set args=-i 1.txt -o out.txt

set CD=%~dp0

rem можно указать, где лежит JDK, если без этого не находится команда java
::set JAVA_HOME=C:\Program Files\Java\jdk-16.0.2\bin

set JAVA=java
if not "%JAVA_HOME%"=="" (
  set JAVA="%JAVA_HOME%\bin\%JAVA%"
)

rem set test=1
rem set wiтdow=--window
rem set io=-i .\input%test%.txt -o .\output%test%.txt


rem В данном случае программа запустится даже если не будет
rem всех перечисленных jar-файлов
rem (т.к. библиотека forms_rt подключалась только для того,
rem  чтобы в режиме генерации кода форм из gui-дизайнера IDEA
rem  не подсвечивала, как несуществующие, некоторые классы,
rem  но при компиляции IDEA код необходимых библиотек включает
rem  в проект)
rem
rem Также обратите внимание на .\out\artifacts\jar
rem

set M2_HOME=%USERPROFILE%\.m2
set CP=%CD%\out\production\task7(13)
set CP=%CP%;%M2_HOME%\repository\com\intellij\forms_rt\7.0.3\forms_rt-7.0.3.jar
set CP=%CP%;%M2_HOME%\repository\asm\asm-commons\3.0\asm-commons-3.0.jar
set CP=%CP%;%M2_HOME%\repository\asm\asm-tree\3.0\asm-tree-3.0.jar
set CP=%CP%;%M2_HOME%\repository\asm\asm\3.0\asm-3.0.jar
set CP=%CP%;%M2_HOME%\repository\com\jgoodies\forms\1.1-preview\forms-1.1-preview.jar
set CP=%CP%;%M2_HOME%\repository\jdom\jdom\1.0\jdom-1.0.jar
set CP=%CP%;%M2_HOME%\repository\commons-cli\commons-cli\1.4\commons-cli-1.4.jar

%JAVA% -classpath "%CP%" ru.vsu.cs.course1.Program %args%