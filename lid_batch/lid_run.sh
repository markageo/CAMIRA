#!/bin/bash

# -------------------------------------------------------------------------------------------------------------- #

# Function to perform division with decimal precision and add a leading zero if the result is less than 1
# Arguments:
#   $1: Numerator
#   $2: Denominator
#   $3: Decimal precision
function divide_with_precision() {
    if [ $# -ne 3 ]; then
        echo "Usage: divide_with_precision <numerator> <denominator> <precision>"
        return 1
    fi

    local numerator="$1"
    local denominator="$2"
    local precision="$3"

    # Use 'scale' in bc to set the decimal precision
    local result=$(echo "scale=$precision; $numerator / $denominator" | bc)

    # Check if the result is less than 1, and if so, add a leading zero
    if (( $(echo "$result < 1" | bc -l) )); then
        result="0$result"
    fi

    echo "$result" 
}

# -------------------------------------------------------------------------------------------------------------- #



inputFileExtension=".inp"
inputFilename="lid"$inputFileExtension
inputFileTemplateFilename="lid_template"$inputFileExtension

# Log file for output from solver
logFileName="solver_terminal.log"
> $logFileName

# Mesh
meshesFolderName="meshes"
# meshes=("22x22x22" 
#         "46x46x46" 
#         "100x100x100"
#         "216x216x216")
meshes=("22x22x22" 
        "46x46x46")
meshFilePrefix="mesh_uniform"

reynoldsNumbers=(200 1000)

# List of sweeping directions
sweepDirections=("+x" "+y" "+z")

# Filename for the rest of the solver properties

# Check for results folder and create it if it doesn't exist
resultsFolderName="results"
if [ ! -d $resultsFolderName ]; then
    mkdir $resultsFolderName
fi

# Loop through Re numbers
for reIdx in ${!reynoldsNumbers[@]}; do

    # Loop through each mesh
    for meshIdx in ${!meshes[@]}; do

        # Plane sweep directions
        for planeSweepDirIdx in ${!sweepDirections[@]}; do

                # Line sweep directions
                for lineSweepDirIdx in ${!sweepDirections[@]}; do

                    re=${reynoldsNumbers[$reIdx]}
                    nu=$(divide_with_precision 1 $re 8)
                    mesh=${meshes[$meshIdx]}
                    planeSweepDir=${sweepDirections[$planeSweepDirIdx]}
                    lineSweepDir=${sweepDirections[$lineSweepDirIdx]}

                    # Line sweep direction cannot be same as plane sweep direction
                    if [ $lineSweepDir == $planeSweepDir ]; then
                        continue
                    fi

                    # Check if the output directory exists, make sure it is empty
                    caseDirectoryName=$resultsFolderName/"lid_re"$re"_"$mesh"_p"$planeSweepDir"_l"$lineSweepDir

                    if [ -d $caseDirectoryName ]; then
                        rm -r $caseDirectoryName/* 2> /dev/null
                    else
                        mkdir $caseDirectoryName
                    fi

                    # Create the input file and put it in the output directory 
                    meshFilename=$meshesFolderName/$meshFilePrefix"_"$mesh$inputFileExtension

                    # Output file names
                    residualHistoryFilename=$caseDirectoryName/residualHistory.csv
                    profilingFilename=$caseDirectoryName/profiling.dat
                    fieldOutputFilename=$caseDirectoryName/fields.vtk


                    # Create and fill input file from template
                    cp $inputFileTemplateFilename $inputFilename

                    # Viscosity
                    sed -i -e 's!$NU!'"$nu"'!g' $inputFilename

                    # Mesh file
                    sed -i -e 's!$MESH_FILENAME!'"$meshFilename"'!g' $inputFilename

                    # Sweep directions
                    sed -i -e 's!$PLANE_SWEEP_DIRECTION!'"$planeSweepDir"'!g' $inputFilename
                    sed -i -e 's!$LINE_SWEEP_DIRECTION!'"$lineSweepDir"'!g' $inputFilename

                    # Output files
                    sed -i -e 's!$RESIDUAL_HISTORY_FILENAME!'"$residualHistoryFilename"'!g' $inputFilename
                    sed -i -e 's!$PROFILING_FILENAME!'"$profilingFilename"'!g' $inputFilename
                    sed -i -e 's!$FIELD_OUTPUT_FILENAME!'"$fieldOutputFilename"'!g' $inputFilename

                    # Copy the input file into the case directory for records
                    cp $inputFilename $caseDirectoryName/$inputFilename


                    # Header 
                    echo $'\n\n' >> $logFileName
                    echo "================= Now running case ================= " | tee -a $logFileName
                    echo "Re:                    "$re             | tee -a $logFileName
                    echo "Mesh:                  "$mesh           | tee -a $logFileName
                    echo "Plane sweep direction: "$planeSweepDir  | tee -a $logFileName
                    echo "Line sweep direction:  "$lineSweepDir   | tee -a $logFileName
                    
                    # Run the code with the input file, need to press enter after meshing has been complete to start solver
                    yes '' | ./main $inputFilename >> $logFileName

                    echo "               * * * Complete! * * *                 " 
                    echo -e

            done

        done

    done

done

rm $inputFilename
