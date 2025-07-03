import './index.css';
import FastaUploader from "./components/Fastauploader";
import AnalysisPage from "./components/AnalysisPage";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<FastaUploader />} />
        <Route path="/analysis" element={<AnalysisPage />} />
      </Routes>
    </Router>
  );
}

export default App;