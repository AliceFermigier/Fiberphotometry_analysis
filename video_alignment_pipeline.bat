@echo off

cd /d "%~dp0"

echo Activating virtual environment
call ".venv\Scripts\activate.bat"

echo Running video_alignment.py
python "modules\behaviour\video_alignment.py"

pause