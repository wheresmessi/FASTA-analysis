from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
from validator import validate_fasta_file
from werkzeug.utils import secure_filename
from Bio import SeqIO
import os

app = Flask(__name__, static_folder="../frontend/build", static_url_path="/")
CORS(app)

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

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

    return jsonify({
        "valid_count": len(valid_records),
        "invalid_count": len(invalids),
        "invalid_details": invalids,
        "valid_file_path": valid_file_path if valid_records else None
    })

# Serve React frontend
@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve_react(path):
    if path != "" and os.path.exists(app.static_folder + "/" + path):
        return send_from_directory(app.static_folder, path)
    else:
        return send_from_directory(app.static_folder, "index.html")

if __name__ == "__main__":
    app.run(debug=True)