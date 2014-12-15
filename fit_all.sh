#bash

for D in 'find . -type d'
do
    python -m analysis_scripts.Jac_test_module
done
