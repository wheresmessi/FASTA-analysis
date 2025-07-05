import React, { useEffect, useState } from "react";
import { useLocation } from "react-router-dom";
import Papa from "papaparse";

function ResultsPage() {
  const location = useLocation();
  const { csv_url, txt_url, message } = location.state || {};
  const [tableData, setTableData] = useState([]);
  const [headers, setHeaders] = useState([]);

  // Fetch and parse CSV
  useEffect(() => {
    if (csv_url) {
      fetch(`http://localhost:5000${csv_url}`)
        .then((res) => res.text())
        .then((csv) => {
          const parsed = Papa.parse(csv, { header: true });
          setTableData(parsed.data);
          setHeaders(parsed.meta.fields || []);
        })
        .catch((err) => {
          console.error("CSV parsing failed:", err);
        });
    }
  }, [csv_url]);

  return (
    <div className="fasta-container">
      <h2 className="fasta-title">Analysis Results</h2>
      {message && <p className="fasta-valid">{message}</p>}

      <div style={{ marginTop: "24px" }}>
        <a
          href={`http://localhost:5000${csv_url}`}
          target="_blank"
          rel="noopener noreferrer"
          className="fasta-upload-btn"
        >
          View CSV Summary
        </a>
        
        <a
          href={`http://localhost:5000${txt_url}`}
          target="_blank"
          rel="noopener noreferrer"
          className="fasta-upload-btn"
          style={{ marginTop: "16px" }}
        >
          View TXT Summary
        </a>
      </div>

      {tableData.length > 0 && (
        <div style={{ marginTop: "32px", overflowX: "auto" }}>
          <h3 style={{ marginBottom: "12px" }}>Preview Table</h3>
          <table className="table-auto border-collapse border w-full text-sm">
            <thead>
              <tr>
                {headers.map((header, idx) => (
                  <th key={idx} className="border px-3 py-2 bg-gray-100">{header}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {tableData.map((row, rowIdx) => (
                <tr key={rowIdx}>
                  {headers.map((h, colIdx) => (
                    <td key={colIdx} className="border px-3 py-2">{row[h]}</td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

export default ResultsPage;