#!/bin/bash


inputFileExtension=".inp"
inputFilename="CEDVAL-A1-7"$inputFileExtension
inputFileTemplateFilename="CEDVAL-A1-7-template"$inputFileExtension
executableDir="/home/mark/PhD/Code/CAMIRA"
outputDir="/home/mark/PhD/Code/CAMIRA/output"

# Mesh
meshNumbers=(1 2 3 4)
nCoarseLevelsArray=(2 3 3 4)

# List of sweeping directions
sweepDirections=("+x" "+y" "+z")

# Check for results folder and create it if it doesn't exist
resultsFolderName=$outputDir/"CEDVAL-A1-7"
if [ ! -d $resultsFolderName ]; then
    mkdir $resultsFolderName
fi

# Log file for output from solver
logFileName=$resultsFolderName/"solver_terminal.log"
> $logFileName


# Loop through each mesh
for meshIdx in ${!meshNumbers[@]}; do

    # Plane sweep directions
    for planeSweepDirIdx in ${!sweepDirections[@]}; do

            # Line sweep directions
            for lineSweepDirIdx in ${!sweepDirections[@]}; do

                meshNumber=${meshNumbers[$meshIdx]}
                nCoarseLevels=${nCoarseLevelsArray[$meshIdx]}
                planeSweepDir=${sweepDirections[$planeSweepDirIdx]}
                lineSweepDir=${sweepDirections[$lineSweepDirIdx]}

                # Line sweep direction cannot be same as plane sweep direction
                if [ $lineSweepDir == $planeSweepDir ]; then
                    continue
                fi

                # Check if the output directory exists, make sure it is empty
                caseDirectoryName=$resultsFolderName/"mesh_"$meshNumber"_p"$planeSweepDir"_l"$lineSweepDir

                if [ -d $caseDirectoryName ]; then
                    rm -r $caseDirectoryName/* 2> /dev/null
                else
                    mkdir $caseDirectoryName
                fi

                # Create the input file and put it in the output directory 
                meshFilename=$meshesFolderName/$meshFilePrefix"_"$mesh$inputFileExtension

                # Create and fill input file from template
                cp $inputFileTemplateFilename $inputFilename

                # Mesh file
                sed -i -e 's!$MESH_FILE_NUMBER!'"$meshNumber"'!g' $inputFilename

                # Sweep directions
                sed -i -e 's!$PLANE_SWEEP_DIRECTION!'"$planeSweepDir"'!g' $inputFilename
                sed -i -e 's!$LINE_SWEEP_DIRECTION!'"$lineSweepDir"'!g' $inputFilename

                # Number of coarse levels
                sed -i -e 's!$COARSE_LEVELS!'"$nCoarseLevels"'!g' $inputFilename

                # Case directory
                sed -i -e 's!$CASE_DIRECTORY!'"$caseDirectoryName"'!g' $inputFilename

                # Copy the input file into the case directory for records
                cp $inputFilename $caseDirectoryName/$inputFilename


                # Header 
                echo $'\n\n' >> $logFileName
                echo "================= CEDVAL A1-7 ================= " | tee -a $logFileName
                echo "Mesh:                  "$meshNumber               | tee -a $logFileName
                echo "Plane sweep direction: "$planeSweepDir            | tee -a $logFileName
                echo "Line sweep direction:  "$lineSweepDir             | tee -a $logFileName
                
                # Run the code with the input file
                $executableDir/camira $inputFilename >> $logFileName

                echo "               * * * Complete! * * *                 " 
                echo -e

        done

    done

done


rm $inputFilename
