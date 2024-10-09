GROMACS Input Generator with UFF and Mixed Force Fields

This Python script automates the generation of molecular connectivity and UFF (Universal Force Field) parameters, converting them into GROMACS-compatible input files using Open Babel's OBGMX methods. It processes .mdp and .car files generated via Material Studio.

Features:

Automatically generates molecular connectivity and applies the UFF force field.
Supports mixing of force fields (OPLS, CVFF, and UFF).
Converts force field parameters to GROMACS-friendly formats.


How to Use:

Ensure that all required files (.mdp, .car) are placed inside a folder named after the structure you are working with.
Set the force_field_mixing parameter:
If force_field_mixing=ON, the script will apply OPLS and CVFF force fields where applicable, while using UFF for any remaining undefined connections.

