@chcp 65001 
:: Изменение кодировки и прочая фигня, забейте

:: Предварительная настройка путей

:: Путь до Java необходимой версии 
@set javaDir="C:\Program Files\Java\jdk-16.0.2\bin"
@set PATH=%javaDir%;%PATH%
:: Путь до корня проекта (out в нем не просто так)
@set cp=C:\Users\Kate\Documents\tasks\8_19int\Arr2JTableDemo\out\production\Arr2JTableDemo
:: Фигня с первой строки каждого java файла + название "главного" файла
@set main=ru.vsu.cs.course1.Program

:: Готовые аргументы 
::
::
:: МЕНЯТЬ ТУТ АРГУМЕНТЫ 
::
::
@set test=1
@set wiтdow=--window
@set io=-i .\input%test%.txt -o .\output%test%.txt

:: Сама команда
java -classpath %cp% %main% %io%
::java -classpath %cp% %main% %wiтdow%

:: Ввод аргументов с консоли (упрощение и сокращение команды
::set /p args="Enter args: "
::java -classpath %cp% %main% %args%

:: Ввод пустой строки для предотвращения закрытия командной строки после выполнения
set /p apgspgpaslg=""