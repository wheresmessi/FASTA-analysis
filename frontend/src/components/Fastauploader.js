import React, { useState } from "react";
import axios from "axios";
import { useNavigate } from "react-router-dom";

const API_BASE = process.env.REACT_APP_API_URL;

function FastaUploader() {
  const [file, setFile] = useState(null);
  const [result, setResult] = useState(null);
  const navigate = useNavigate();

  const handleUpload = async () => {
    const formData = new FormData();
    formData.append("file", file);
    try {
      const res = await axios.post(`${API_BASE}/upload`, formData);
      setResult(res.data);
    } catch (error) {
      console.error("Upload failed:", error);
    }
  };

  const handleAnalysis = () => {
    navigate("/analysis", { state: { validFilePath: result.valid_file_path } });
  };

  return (
    <div className="fasta-container">
      <h2 className="fasta-title">Upload FASTA File</h2>
      <div className="fasta-form">
        <input type="file" className="fasta-input" onChange={(e) => setFile(e.target.files[0])} />
        <button onClick={handleUpload} className="fasta-upload-btn">Upload</button>
      </div>
      {result && (
        <div className="fasta-results">
          <h3>Results:</h3>
          <p className="fasta-valid">✅ <strong>Valid Records:</strong> {result.valid_count}</p>
          <button className="fasta-upload-btn" style={{ marginTop: "20px" }} onClick={handleAnalysis}>
            Perform Analysis
          </button>
        </div>
      )}
    </div>
  );
}

export default FastaUploader;