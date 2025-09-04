#!/bin/bash


inputFileExtension=".inp"
inputFilename="CEDVAL-B1-1"$inputFileExtension
inputFileTemplateFilename="CEDVAL-B1-1-template"$inputFileExtension
executableDir="/home/mark/PhD/Code/CAMIRA"
outputDir="/home/mark/PhD/Code/CAMIRA/output"


# List of sweeping directions
sweepDirections=("+x" "+y" "+z")

# Check for results folder and create it if it doesn't exist
resultsFolderName=$outputDir/"CEDVAL-B1-1"
if [ ! -d $resultsFolderName ]; then
    mkdir $resultsFolderName
fi

# Log file for output from solver
logFileName=$resultsFolderName/"solver_terminal.log"
> $logFileName


# Plane sweep directions
for planeSweepDirIdx in ${!sweepDirections[@]}; do

        # Line sweep directions
        for lineSweepDirIdx in ${!sweepDirections[@]}; do

            planeSweepDir=${sweepDirections[$planeSweepDirIdx]}
            lineSweepDir=${sweepDirections[$lineSweepDirIdx]}

            # Line sweep direction cannot be same as plane sweep direction
            if [ $lineSweepDir == $planeSweepDir ]; then
                continue
            fi

            # Check if the output directory exists, make sure it is empty
            caseDirectoryName=$resultsFolderName/"p"$planeSweepDir"_l"$lineSweepDir

            if [ -d $caseDirectoryName ]; then
                rm -r $caseDirectoryName/* 2> /dev/null
            else
                mkdir $caseDirectoryName
            fi

            # Create the input file and put it in the output directory 
            meshFilename=$meshesFolderName/$meshFilePrefix"_"$mesh$inputFileExtension

            # Create and fill input file from template
            cp $inputFileTemplateFilename $inputFilename

            # Sweep directions
            sed -i -e 's!$PLANE_SWEEP_DIRECTION!'"$planeSweepDir"'!g' $inputFilename
            sed -i -e 's!$LINE_SWEEP_DIRECTION!'"$lineSweepDir"'!g' $inputFilename

            # Case directory
            sed -i -e 's!$CASE_DIRECTORY!'"$caseDirectoryName"'!g' $inputFilename

            # Copy the input file into the case directory for records
            cp $inputFilename $caseDirectoryName/$inputFilename


            # Header 
            echo $'\n\n' >> $logFileName
            echo "================= CEDVAL B1-1 ================= " | tee -a $logFileName
            echo "Plane sweep direction: "$planeSweepDir            | tee -a $logFileName
            echo "Line sweep direction:  "$lineSweepDir             | tee -a $logFileName
            
            # Run the code with the input file
            $executableDir/camira $inputFilename >> $logFileName

            echo "               * * * Complete! * * *                 " 
            echo -e

    done

done



rm $inputFilename
