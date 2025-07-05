import React from "react";
import { useLocation, useNavigate } from "react-router-dom";
import axios from "axios";

const API_BASE = process.env.REACT_APP_API_URL;

function AnalysisPage() {
  const location = useLocation();
  const navigate = useNavigate();
  const { validFilePath } = location.state || {};

  const handleGCAnalysis = async () => {
    try {
      const response = await axios.post(`${API_BASE}/run-analysis`, { valid_file_path: validFilePath });
      navigate("/results", { state: response.data });
    } catch (err) {
      alert("Analysis failed!");
      console.error(err);
    }
  };

  return (
    <div className="fasta-container">
      <h2 className="fasta-title">Analysis Page</h2>
      <div className="analysis-buttons" style={{ display: "flex", flexDirection: "column", gap: "16px", marginTop: "32px" }}>
        <button className="fasta-upload-btn" onClick={handleGCAnalysis}>GC content & nucleotide content</button>
        <button className="fasta-upload-btn">ENC plot</button>
        <button className="fasta-upload-btn">Correspondence Analysis (COA)</button>
        <button className="fasta-upload-btn">RSCU / optimal codon usage</button>
        <button className="fasta-upload-btn">PR2 bias</button>
        <button className="fasta-upload-btn">Neutrality plot</button>
      </div>
    </div>
  );
}

export default AnalysisPage;