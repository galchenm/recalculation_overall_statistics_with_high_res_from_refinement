#!/usr/bin/env python3

import os
import sys
import subprocess
import time
import argparse
import tkinter as tk
from tkinter import ttk
import pandas as pd

REFRESH_INTERVAL_MS = 3000  # Refresh every 3 seconds

class CSVViewerApp:
    def __init__(self, root, csv_path):
        self.root = root
        self.csv_path = csv_path
        self.dataframe = pd.DataFrame()

        self.setup_ui()
        self.wait_for_file()  # Check for the file before proceeding
        self.schedule_refresh()

    def setup_ui(self):
        self.root.title("Live CSV Viewer")
        self.root.geometry("1000x600")

        # Table
        self.tree = ttk.Treeview(self.root, show="headings")
        self.tree.pack(fill=tk.BOTH, expand=True)

        # Scrollbars
        scrollbar_y = ttk.Scrollbar(self.tree, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar_y.set)
        scrollbar_y.pack(side="right", fill="y")

        scrollbar_x = ttk.Scrollbar(self.tree, orient="horizontal", command=self.tree.xview)
        self.tree.configure(xscroll=scrollbar_x.set)
        scrollbar_x.pack(side="bottom", fill="x")

        # Style
        style = ttk.Style()
        style.configure("Treeview", rowheight=25)
        style.map("Treeview", background=[('selected', '#ececec')])
        style.configure("Treeview.Heading", font=('Helvetica', 10, 'bold'))

    def wait_for_file(self, timeout_hours=6, check_interval=2000):
        # Wait until the file appears
        timeout_seconds = timeout_hours * 3600
        start_time = time.time()
        
        while not os.path.isfile(self.csv_path):
            elapsed = time.time() - start_time
            if elapsed > timeout_seconds:
                print(f"[ERROR] File not found at {self.csv_path} after {timeout_hours} hours. Exiting.")
                sys.exit(1)
            print(f"[INFO] Waiting for file: {self.csv_path}")
            time.sleep(check_interval)  # Wait for check_interval seconds

        print(f"[INFO] File found: {self.csv_path}")
        self.load_and_display_csv()  # Load the file once it is found

    def load_and_display_csv(self):
        # Try to load the CSV file
        try:
            df = pd.read_csv(self.csv_path, delimiter=';')
            self.dataframe = df
        except Exception as e:
            print(f"[ERROR] Failed to load CSV: {e}")
            self.dataframe = pd.DataFrame()  # Empty DataFrame on error

        self.update_table()

    def update_table(self):
        # Update the table with the new DataFrame
        df_to_display = self.dataframe

        self.tree.delete(*self.tree.get_children())  # Clear existing rows
        self.tree["columns"] = []  # Reset the columns

        if df_to_display.empty:
            return

        # Set new columns based on the DataFrame
        columns = list(df_to_display.columns)
        self.tree["columns"] = columns

        for col in columns:
            self.tree.heading(col, text=col)  # Set column headers
            self.tree.column(col, width=250, anchor="w")  # Set column widths

        # Insert new rows
        for _, row in df_to_display.iterrows():
            self.tree.insert("", "end", values=list(row))

    def schedule_refresh(self):
        # Refresh data every REFRESH_INTERVAL_MS milliseconds
        self.load_and_display_csv()
        self.root.after(REFRESH_INTERVAL_MS, self.schedule_refresh)


def main():
    parser = argparse.ArgumentParser(description="Live CSV Viewer")
    parser.add_argument("--csv", required=True, help="Path to CSV file")
    args = parser.parse_args()

    # Wait until the CSV file exists before proceeding
    while not os.path.isfile(args.csv):
        print(f"[ERROR] File '{args.csv}' not found.")
        time.sleep(5)  # Wait for the file to be created

    # Initialize the application
    root = tk.Tk()
    app = CSVViewerApp(root, args.csv)
    root.mainloop()


if __name__ == "__main__":
    main()
