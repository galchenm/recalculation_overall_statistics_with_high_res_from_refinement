#!/usr/bin/env python3

import subprocess
import sys
import time
import io
import os

# === Ensure Flask is installed ===
try:
    import flask
except ImportError:
    print("[INFO] Flask not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "flask"])
    import flask  # retry import after installation

# === Ensure portalocker is installed for file locking ===
try:
    import portalocker
except ImportError:
    print("[INFO] portalocker not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "portalocker"])
    import portalocker

from flask import Flask, render_template_string, jsonify
import argparse
import csv

app = Flask(__name__)
CSV_PATH = ""

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Live CSV Viewer</title>
    <style>
        body { font-family: sans-serif; padding: 1em; background: #f9f9f9; }
        table { border-collapse: collapse; width: 100%; margin-top: 10px; }
        th, td { border: 1px solid #ccc; padding: 6px 10px; text-align: left; }
        th { background-color: #eee; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        #table-container { height: 70vh; overflow-y: auto; border: 1px solid #ccc; }
    </style>
</head>
<body>
    <h1>Live CSV Viewer</h1>
    <div id="table-container">
        <table id="csv-table">
            <thead></thead>
            <tbody></tbody>
        </table>
    </div>

    <script>
        function fetchCSV() {
            fetch('/data')
                .then(response => response.json())
                .then(data => {
                    const thead = document.querySelector("#csv-table thead");
                    const tbody = document.querySelector("#csv-table tbody");

                    // Clear current content
                    thead.innerHTML = "";
                    tbody.innerHTML = "";

                    if (data.length === 0) return;

                    // Header
                    const headerRow = document.createElement("tr");
                    Object.keys(data[0]).forEach(key => {
                        const th = document.createElement("th");
                        th.textContent = key;
                        headerRow.appendChild(th);
                    });
                    thead.appendChild(headerRow);

                    // Rows
                    data.forEach(row => {
                        const tr = document.createElement("tr");
                        Object.values(row).forEach(val => {
                            const td = document.createElement("td");
                            td.textContent = val;
                            tr.appendChild(td);
                        });
                        tbody.appendChild(tr);
                    });

                    autoScroll();
                });
        }

        function autoScroll() {
            const container = document.getElementById("table-container");
            container.scrollTop = container.scrollHeight;
        }

        setInterval(fetchCSV, 3000);  // Refresh every 3 seconds
        window.onload = fetchCSV;
    </script>
</body>
</html>
"""

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/data')
def data():
    rows = []
    max_attempts = 3
    attempt = 0

    while attempt < max_attempts:
        try:
            if os.path.isfile(CSV_PATH):
                with open(CSV_PATH, "r", newline='') as f:
                    portalocker.lock(f, portalocker.LOCK_SH)  # Shared lock for reading
                    content = f.read()
                    portalocker.unlock(f)
                    reader = csv.DictReader(io.StringIO(content), delimiter=';')
                    for row in reader:
                        clean_row = {k: (v if v is not None else "") for k, v in row.items()}
                        rows.append(clean_row)
            break
        except (csv.Error, OSError) as e:
            attempt += 1
            time.sleep(0.2)  # brief pause before retry
            if attempt == max_attempts:
                print(f"[WARNING] Failed to read CSV after {max_attempts} attempts: {e}")
    return jsonify(rows)

def wait_for_file(csv_path, timeout_hours=6, check_interval=2):
    timeout_seconds = timeout_hours * 3600
    start_time = time.time()

    print(f"[INFO] Waiting for file: {csv_path}")

    while not os.path.isfile(csv_path):
        elapsed = time.time() - start_time
        if elapsed > timeout_seconds:
            print(f"[ERROR] File not found at {csv_path} after {timeout_hours} hours. Exiting.")
            sys.exit(1)
        time.sleep(check_interval)  # check every 30 seconds

    print(f"[INFO] File found: {csv_path}")

def main():
    global CSV_PATH
    parser = argparse.ArgumentParser(description="Live CSV Table Viewer")
    parser.add_argument('--csv', required=True, help='Path to the CSV file to monitor')
    parser.add_argument('--port', type=int, default=5000, help='Port to serve on (default: 5000)')
    args = parser.parse_args()

    CSV_PATH = args.csv
    print(157)
    wait_for_file(csv_path=CSV_PATH)

    print(f"Serving live view of '{CSV_PATH}' at http://localhost:{args.port}")
    app.run(debug=False, use_reloader=False, port=args.port)

if __name__ == '__main__':
    main()
