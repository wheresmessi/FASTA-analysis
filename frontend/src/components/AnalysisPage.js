import React from "react";

function AnalysisPage() {
  return (
    <div className="fasta-container">
      <h2 className="fasta-title">Analysis Page</h2>
      <div className="analysis-buttons" style={{ display: "flex", flexDirection: "column", gap: "16px", marginTop: "32px" }}>
        <button className="fasta-upload-btn">GC content &amp; nucleotide content</button>
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