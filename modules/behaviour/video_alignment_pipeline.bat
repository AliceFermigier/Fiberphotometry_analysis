@echo off
echo === Activating venv ===
call c:\Users\alice\Documents\GitHub\Fiberphotometry_analysis\.venv\Scripts\activate.bat

REM Going to project root
cd /d c:\Users\alice\Documents\GitHub\Fiberphotometry_analysis

echo === Running video_alignment.py ===
python modules\behaviour\video_alignment.py

pause