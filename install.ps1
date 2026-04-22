$ErrorActionPreference = "Stop"
if ($env:PYTHON_BIN) {
    & $env:PYTHON_BIN install.py @args
} else {
    py -3 install.py @args
}
