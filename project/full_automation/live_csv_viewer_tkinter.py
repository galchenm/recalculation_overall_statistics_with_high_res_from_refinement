#!/usr/bin/env python3

import os
import sys
import subprocess

# Ensure pandas is installed
try:
    import pandas as pd
except ImportError:
    print("[INFO] pandas not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd

import argparse
import tkinter as tk
from tkinter import ttk

REFRESH_INTERVAL_MS = 3000  # Refresh every 3 seconds

class CSVViewerApp:
    def __init__(self, root, csv_path):
        self.root = root
        self.csv_path = csv_path
        self.search_var = tk.StringVar()

        self.root.title("Live CSV Viewer")
        self.root.geometry("1000x600")

        # Search field
        search_frame = tk.Frame(root)
        search_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Label(search_frame, text="Search:").pack(side=tk.LEFT)
        tk.Entry(search_frame, textvariable=self.search_var).pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.search_var.trace_add("write", lambda *_: self.update_table())

        # Treeview (table)
        self.tree = ttk.Treeview(root, show="headings")
        self.tree.pack(fill=tk.BOTH, expand=True)

        # Scrollbar
        scrollbar_y = ttk.Scrollbar(self.tree, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar_y.set)
        scrollbar_y.pack(side="right", fill="y")

        self.dataframe = pd.DataFrame()
        self.load_and_display_csv()
        self.schedule_refresh()

    def load_and_display_csv(self):
        if not os.path.exists(self.csv_path):
            return

        try:
            df = pd.read_csv(self.csv_path)
            self.dataframe = df
        except Exception as e:
            print(f"[ERROR] Failed to load CSV: {e}")
            return

        self.update_table()

    def update_table(self):
        search_term = self.search_var.get().lower()
        df = self.dataframe

        if search_term:
            df = df[df.apply(lambda row: row.astype(str).str.lower().str.contains(search_term).any(), axis=1)]

        # Clear table
        for col in self.tree["columns"]:
            self.tree.heading(col, text="")
        self.tree.delete(*self.tree.get_children())

        # Set new columns
        if not df.empty:
            self.tree["columns"] = list(df.columns)
            for col in df.columns:
                self.tree.heading(col, text=col)

            for _, row in df.iterrows():
                self.tree.insert("", "end", values=list(row))

    def schedule_refresh(self):
        self.load_and_display_csv()
        self.root.after(REFRESH_INTERVAL_MS, self.schedule_refresh)


def main():
    parser = argparse.ArgumentParser(description="Live CSV Viewer (tkinter)")
    parser.add_argument("--csv", required=True, help="Path to CSV file")
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        print(f"[ERROR] File '{args.csv}' does not exist.")
        sys.exit(1)

    root = tk.Tk()
    app = CSVViewerApp(root, args.csv)
    root.mainloop()

if __name__ == "__main__":
    main()
