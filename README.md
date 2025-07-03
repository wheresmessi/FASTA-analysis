# FASTA Uploader & Analysis Web App

This project is a full-stack web application for uploading FASTA files, validating coding sequences, and performing various sequence analyses.

## Features

- Upload and validate FASTA files (checks for valid CDS: length, start/stop codons, etc.)
- Displays count of valid records
- "Perform Analysis" button redirects to an analysis page
- Analysis page with buttons for:
  - GC content & nucleotide content
  - ENC plot
  - Correspondence Analysis (COA)
  - RSCU / optimal codon usage
  - PR2 bias
  - Neutrality plot

## Tech Stack

- **Frontend:** React (with React Router), Axios, CSS
- **Backend:** Flask, Biopython

## How to Run

### 1. Clone the repository

```sh
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```

### 2. Backend Setup

```sh
cd backend
pip install -r requirements.txt
# If requirements.txt doesn't exist, install manually:
pip install flask flask-cors biopython
python app.py
```

### 3. Frontend Setup

Open a new terminal and run:

```sh
cd frontend
npm install
npm start
```

- The React app runs on [http://localhost:3000](http://localhost:3000)
- The Flask backend runs on [http://localhost:5000](http://localhost:5000)

### 4. Production Build

To serve the React app with Flask, build the frontend:

```sh
cd frontend
npm run build
```
Then restart the Flask backend.

## Project Structure

```
backend/
  app.py
  validator.py
frontend/
  src/
    components/
      Fastauploader.js
    AnalysisPage.js
    App.js
    index.css
  package.json
```

## License

MIT

---

**Replace the repo URL with your actual GitHub repository.**  
Add this `README.md` to your project root and commit it:

```sh
git add README.md
git commit -m "Add project README"
git push
```