Get-ChildItem "A:\projects\kansasii\input\*.CSV" | ForEach-Object {
    $tmp = $_.FullName + ".tmp"
    Rename-Item $_.FullName $tmp
    Rename-Item $tmp ($_.BaseName + ".csv")
}





powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

uv python list --all-platforms

uv venv .venv --python 3.14
uv pip install --python .venv/Scripts/python.exe pandas numpy scipy matplotlib pytest -e .

scp -r "A:\projects\kansasii\input\*" mimeul@cluster.s3it.uzh.ch:/shares/sander.imm.uzh/MM/kansasii/data/imm/




.venv\Scripts\python.exe scripts\run_analysis.py `
  --mic-csv "A:\projects\kansasii\input\mic.csv" `
  --screening-map "A:\projects\kansasii\input\screening_map.csv" `
  --out-dir "A:\projects\kansasii\downloads\mic"


.venv\Scripts\Activate.ps1
python scripts\run_analysis.py `
  --mic-csv "A:\projects\kansasii\input\mic.csv" `
  --screening-map "A:\projects\kansasii\input\screening_map.csv" `
  --out-dir "A:\projects\kansasii\downloads\mic"

