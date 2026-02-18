@echo off
cd /d "%~dp0"
echo.
echo Starting Biomarker Pipeline UI...
echo.
echo When the server is ready, open this link in your browser:
echo   http://localhost:8501
echo.
python -m streamlit run app.py --server.headless true
pause
