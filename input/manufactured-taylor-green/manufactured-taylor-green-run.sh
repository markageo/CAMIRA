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
inputFilename="manufactured-taylor-green"$inputFileExtension
inputFileTemplateFilename="manufactured-taylor-green-template"$inputFileExtension

# Log file for output from solver
logFileName="solver_terminal.log"
> $logFileName

# Mesh and advection schemes to use
nCells=(20 40 80 160 320)
advectionSchemes=("upwind" "QUICK")
nCoarseLevels=(0 0 1 2 3)

# Mesh name strings
meshes=()
for m in ${!nCells[@]}; do
    meshes+=(${nCells[$m]}x${nCells[$m]}x${nCells[$m]})
done

# Check for results folder and create it if it doesn't exist
resultsFolderName="results"
if [ ! -d $resultsFolderName ]; then
    mkdir $resultsFolderName
fi

# Loop through advection schemes
for advIdx in ${!advectionSchemes[@]}; do
    advectionScheme=${advectionSchemes[$advIdx]}

    # Loop thorugh meshes
    for meshIdx in ${!meshes[@]}; do
        mesh=${meshes[$meshIdx]}
        
         # Check if the output directory exists, make sure it is empty
        caseDirectoryName=$resultsFolderName/"manufactured_taylor_green_"$mesh"_"$advectionScheme

        if [ -d $caseDirectoryName ]; then
            rm -r $caseDirectoryName/* 2> /dev/null
        else
            mkdir $caseDirectoryName
        fi

        # Create and fill input file from template
        cp $inputFileTemplateFilename $inputFilename

        # Mesh
        sed -i -e 's!$NCELLS!'"${nCells[$meshIdx]}"'!g' $inputFilename

        # Advection scheme
        sed -i -e 's!$ADVECTION_SCHEME!'"$advectionScheme"'!g' $inputFilename

        # Coarse levels
        sed -i -e 's!$COARSE_LEVELS!'"${nCoarseLevels[$meshIdx]}"'!g' $inputFilename

        # Case directories
        sed -i -e 's!$CASE_DIRECTORY!'"$caseDirectoryName"'!g' $inputFilename    


        # Copy the input file into the case directory 
        cp $inputFilename $caseDirectoryName/$inputFilename

        # Header 
        echo $'\n\n' >> $logFileName
        echo "================= Now running case ================= " | tee -a $logFileName
        echo "Advection Scheme:      "$advectionScheme  | tee -a $logFileName
        echo "Mesh:                  "$mesh             | tee -a $logFileName
        
        # Run the code with the input file, need to press enter after meshing has been complete to start solver
        yes '' | ./camira $inputFilename >> $logFileName

        echo "               * * * Complete! * * *                 " 
        echo -e

    done 

done


# rm $inputFilename
