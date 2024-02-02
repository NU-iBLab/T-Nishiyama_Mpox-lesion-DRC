# T-Nishiyama_Mpox-lesion-DRC
This code is an analysis code for the paper "Modeling lesion transition dynamics to clinically characterize mpox patients in the Democratic Republic of the Congo"

## Desciption
Mathematical models describing the lesion and virus　dynamics are stored in monolix files.
See (https://lixoft.com/products/monolix/) for parameter estimation.

The complete dataset used in this study are available from the corresponding author upon reasonable request.

## monolix
'Lesion_dynamics/monolix_4ph.csv' is lesion data of Mpox patient.
'Lesion_dynamics/cure.txt' is the description for mathematical model of lesion dynamics.

'Virus_dynamics/formonolix.csv' is viral load data of Mpox patient.
'Virus_dynamics/simplemodel.txt' is the description for mathematical model of virus dynamics.

These files are used to estimate parameters of the lesion and virus dynamics model in MONOLIX2021R2.
