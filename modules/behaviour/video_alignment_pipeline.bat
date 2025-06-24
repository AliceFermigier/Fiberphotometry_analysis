@echo off
echo Activating virtual environment
call c:\Users\afermigier\Documents\GitHub\Fiberphotometry_analysis\.venv\Scripts\activate.bat

REM Going to project root
cd /d c:\Users\afermigier\Documents\GitHub\Fiberphotometry_analysis

echo Running video_alignment.py
python modules\behaviour\video_alignment.py

pause