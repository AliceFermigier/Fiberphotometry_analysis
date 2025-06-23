@echo off
echo === Activating venv ===
call c:\Users\afermigier\Documents\GitHub\Fiberphotometry_analysis\.venv\Scripts\activate.bat

echo === Go to the project root ===
cd /d c:\Users\afermigier\Documents\GitHub\Fiberphotometry_analysis

echo === Running loader.py ===
python scripts\loader.py

echo === Running video_alignment.py ===
python modules\behaviour\video_alignment.py

pause