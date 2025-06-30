import pytest
import tempfile
import csv
import time
import requests
import os
import sys
import multiprocessing
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from life_csv_web_viewer import app, wait_for_file  # replace 'your_script' with your actual filename (no .py)

# Enable test mode
app.config["TESTING"] = True

# --- Helper to run Flask app in a separate process ---
def run_app(csv_path, port):
    from life_csv_web_viewer import CSV_PATH
    CSV_PATH = csv_path  # overwrite global path
    wait_for_file(CSV_PATH, timeout_hours=0.001, check_interval=0.1)
    app.run(port=port, debug=False, use_reloader=False)

# Fixture to create a temporary CSV file
@pytest.fixture
def temp_csv_file():
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".csv", newline='', encoding='utf-8') as tmp:
        writer = csv.writer(tmp, delimiter=';')
        writer.writerow(['name', 'value'])  # CSV header
        tmp.flush()
        yield tmp.name
    os.remove(tmp.name)

def test_flask_csv_viewer(temp_csv_file):
    port = 5055  # pick any available port

    # Start the Flask app in a separate process
    proc = multiprocessing.Process(target=run_app, args=(temp_csv_file, port))
    proc.start()
    time.sleep(1.5)  # allow the server to start

    try:
        # Test main page
        index_resp = requests.get(f'http://localhost:{port}/')
        assert index_resp.status_code == 200
        assert "<title>Live CSV Viewer</title>" in index_resp.text

        # Test /data initially returns empty
        data_resp = requests.get(f'http://localhost:{port}/data')
        assert data_resp.status_code == 200
        assert data_resp.json() == []

        # Simulate writing a new row to the CSV
        with open(temp_csv_file, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter=';')
            writer.writerow(['TestName', '123'])
        time.sleep(1.0)  # wait for browser refresh interval

        # Test that /data returns the new row
        data_resp = requests.get(f'http://localhost:{port}/data')
        assert data_resp.status_code == 200
        data = data_resp.json()
        assert len(data) == 1
        assert data[0]['name'] == 'TestName'
        assert data[0]['value'] == '123'

    finally:
        proc.terminate()
        proc.join()
