#!/usr/bin/env python3

import os
import sys
import subprocess
import time
import argparse
import tkinter as tk
from tkinter import ttk

# Ensure pandas is installed
try:
    import pandas as pd
except ImportError:
    print("[INFO] pandas not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd

REFRESH_INTERVAL_MS = 3000000  # Refresh every 3 seconds


class CSVViewerApp:
    def __init__(self, root, csv_path):
        self.root = root
        self.csv_path = csv_path
        self.dataframe = pd.DataFrame()

        self.setup_ui()
        self.load_and_display_csv()
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

        # Style
        style = ttk.Style()
        style.configure("Treeview", rowheight=25)
        style.map("Treeview", background=[('selected', '#ececec')])
        style.configure("Treeview.Heading", font=('Helvetica', 10, 'bold'))

    def load_and_display_csv(self):
        if not os.path.exists(self.csv_path):
            print(f"[WARN] File not found: {self.csv_path}")
            return

        for attempt in range(3):
            try:
                df = pd.read_csv(self.csv_path, delimiter=';')
                self.dataframe = df
                break
            except Exception as e:
                print(f"[ERROR] Failed to load CSV (attempt {attempt + 1}): {e}")
                time.sleep(3000)
        else:
            print("[ERROR] Could not read CSV after multiple attempts.")
            self.dataframe = pd.DataFrame()

        self.update_table()

    def update_table(self):
        df_to_display = self.dataframe

        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = []

        if df_to_display.empty:
            return

        columns = list(df_to_display.columns)
        self.tree["columns"] = columns

        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=250, anchor="w")

        for _, row in df_to_display.iterrows():
            self.tree.insert("", "end", values=list(row))

    def schedule_refresh(self):
        self.load_and_display_csv()
        self.root.after(REFRESH_INTERVAL_MS, self.schedule_refresh)


def main():
    parser = argparse.ArgumentParser(description="Live CSV Viewer")
    parser.add_argument("--csv", required=True, help="Path to CSV file")
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        print(f"[ERROR] File '{args.csv}' not found.")
        sys.exit(1)

    root = tk.Tk()
    app = CSVViewerApp(root, args.csv)
    root.mainloop()


if __name__ == "__main__":
    main()
