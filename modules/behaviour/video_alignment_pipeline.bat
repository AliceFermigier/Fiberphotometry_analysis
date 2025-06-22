@echo off
REM Activate your virtual environment
call C:\Users\alice\Documents\GitHub\Fiberphotometry_analysis\.venv\Scripts\activate.bat

REM Go to the project root
cd /d C:\Users\alice\Documents\GitHub\Fiberphotometry_analysis

REM (Optional) Run loader.py if needed
python scripts\loader.py

REM Run your target script
python modules\behaviour\video_alignment.py

pause