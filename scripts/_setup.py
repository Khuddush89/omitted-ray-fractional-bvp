"""Shared paths and helpers for the scripts."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)


def save_json(name, obj):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=False)
    return path


def load_json(name):
    with open(os.path.join(RESULTS, name)) as fh:
        return json.load(fh)


def save_text(name, text):
    path = os.path.join(RESULTS, name)
    with open(path, "w") as fh:
        fh.write(text)
    return path
