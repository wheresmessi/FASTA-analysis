import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import FastaUploader from "./components/Fastauploader";
import AnalysisPage from "./components/AnalysisPage";
import ResultsPage from "./components/ResultsPage";

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<FastaUploader />} />
        <Route path="/analysis" element={<AnalysisPage />} />
        <Route path="/results" element={<ResultsPage />} />
      </Routes>
    </Router>
  );
}

export default App;