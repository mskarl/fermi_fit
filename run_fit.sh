# example bash skript:

source myPythonEnv/bin/activate
conda activate fermi

echo job is running: $1
python run_fit.py --source $1 --force_power_law

