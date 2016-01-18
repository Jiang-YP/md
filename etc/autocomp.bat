@echo off
rem This is a batch processing file to compile all the needed source files into object files.
rem
echo.
echo 1 - Copying the sources from \dev to \src
md ..\tmp
xcopy ..\dev\*.f* ..\tmp\ /y /q
xcopy ..\src\*.f* ..\tmp\ /y /q
cd ..\tmp
echo.
echo 2 - Beginning to compile files
echo.
echo Processing commod.f        ... [01/10]: Common module
call aspcomp commod.f
echo Processing mdvar.f         ... [02/10]: AspenPlus configuration variables definition
call aspcomp mdvar.f
echo Processing d1mach.f        ... [03/10]: Subroutine for machine number, invoked by integal solver
call aspcomp d1mach.f
call aspcomp i1mach.f 
echo Processing dqk15.f         ... [04/10]: Integal solver
call aspcomp dqk15.f
echo Processing diff.f          ... [05/10]: Deriviation solver
call aspcomp diff.f
echo Processing property.f90    ... [06/10]: Subroutines of physical properties
call aspcomp property.f90
echo Processing localflux.f90   ... [07/10]: Subroutine for local heat and mass transfer rates
call aspcomp localflux.f90
echo Processing toolkits.f90    ... [08/10]: toolkits subroutine
call aspcomp toolkits.f90
echo Processing CalcProf.f90    ... [09/10]: Subroutine for get the temperature profiles
call aspcomp CalcProf.f90
echo Processing CalcSOUT.f90    ... [10/10]: Main program
call aspcomp CalcSOUT.f90
echo Process dcmd.f             ... [11/10]: Aspen interface program
call aspcomp dcmd.f
echo.
echo Compilation finished ...
echo.
echo 3 - Linking the object files into dynamic linked library
echo.
call asplink dcmd
echo.
echo 4 - Copying dcmd.dll to the binary directory
xcopy ..\tmp\dcmd.dll ..\bin\ /y /q
xcopy ..\tmp\dcmd.dll "D:\Works\AspenTech\Working Folders\Aspen Plus V7.3\20160105 md with Kim2014 model\" /y /q
type ..\tmp\dcmd.ld
cd ..\etc
echo.
echo 5 - Finished!
