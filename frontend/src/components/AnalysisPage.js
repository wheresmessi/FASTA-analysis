import React from "react";
import { useLocation, useNavigate } from "react-router-dom";
import axios from "axios";

function AnalysisPage() {
  const location = useLocation();
  const navigate = useNavigate();
  const { validFilePath } = location.state || {};

  const handleGCAnalysis = async () => {
    if (!validFilePath) {
      alert("No valid FASTA file path found!");
      return;
    }

    try {
      const response = await axios.post("http://localhost:5000/run-analysis", {
        valid_file_path: validFilePath,
      });

      navigate("/results", { state: response.data });
    } catch (err) {
      console.error("Analysis failed:", err);
      alert("Analysis failed!");
    }
  };

  return (
    <div className="fasta-container">
      <h2 className="fasta-title">Analysis Page</h2>
      <div className="analysis-buttons" style={{ display: "flex", flexDirection: "column", gap: "16px", marginTop: "32px" }}>
        <button className="fasta-upload-btn" onClick={handleGCAnalysis}>
          GC content &amp; nucleotide content
        </button>
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