from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from werkzeug.utils import secure_filename
from Bio import SeqIO
import os
import shutil
from validator import validate_fasta_file
from gc_analysis import run_combined_analysis

# === Absolute Paths ===
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
UPLOAD_FOLDER = os.path.join(BASE_DIR, "uploads")
RESULTS_FOLDER = os.path.join(BASE_DIR, "results")
INPUT_FILE = os.path.join(UPLOAD_FOLDER, "input_for_analysis.txt")
FRONTEND_BUILD = os.path.join(BASE_DIR, "frontend", "build")

# === Flask Setup ===
app = Flask(__name__, static_folder=FRONTEND_BUILD, static_url_path="/")

# ✅ CORS: Allow frontend hosted on Render
CORS(app, resources={r"/*": {"origins": "https://fasta-analysis.onrender.com"}})

# Ensure folders exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

# === Upload and Validate FASTA ===
@app.route("/upload", methods=["POST"])
def upload():
    if 'file' not in request.files:
        return jsonify({"error": "No file uploaded"}), 400

    file = request.files['file']
    filename = secure_filename(file.filename)
    filepath = os.path.join(UPLOAD_FOLDER, filename)
    file.save(filepath)

    valid_records, invalids = validate_fasta_file(filepath)
    valid_file_path = os.path.join(UPLOAD_FOLDER, f"{filename}_valid.fasta")

    if valid_records:
        with open(valid_file_path, "w") as handle:
            SeqIO.write(valid_records, handle, "fasta")
        shutil.copy(valid_file_path, INPUT_FILE)

    return jsonify({
        "valid_count": len(valid_records),
        "invalid_count": len(invalids),
        "invalid_details": invalids,
        "valid_file_path": valid_file_path if valid_records else None
    })

# === Run Analysis ===
@app.route("/run-analysis", methods=["POST"])
def run_analysis():
    if not os.path.exists(INPUT_FILE):
        return jsonify({"error": "Valid FASTA input not found."}), 400

    try:
        run_combined_analysis()
    except Exception as e:
        return jsonify({"error": "Analysis failed", "details": str(e)}), 500

    return jsonify({
        "message": "Analysis complete",
        "csv_url": "/download/combined_summary.csv",
        "txt_url": "/download/combined_summary.txt"
    })

# === Serve Result Files ===
@app.route("/download/<path:filename>", methods=["GET"])
def download_file(filename):
    return send_from_directory(RESULTS_FOLDER, filename, as_attachment=False)

# === Serve React Frontend ===
@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve_react(path):
    target = os.path.join(FRONTEND_BUILD, path)
    if path != "" and os.path.exists(target):
        return send_from_directory(FRONTEND_BUILD, path)
    return send_from_directory(FRONTEND_BUILD, "index.html")

# === Run Flask ===
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)