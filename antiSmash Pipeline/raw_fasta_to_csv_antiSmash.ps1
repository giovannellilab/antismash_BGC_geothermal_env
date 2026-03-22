# Get the directory where the script is being run
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition

# Define the base directory containing all subfolders with .fa files
# This is where you should have your input files
$baseDir = Join-Path -Path $scriptDir -ChildPath "Input"

# Define the main output directory for results (relative to script location)
$outputDir = Join-Path -Path $scriptDir -ChildPath "Pipeline Output"

# Create the main output directory if it does not exist
if (-Not (Test-Path -Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir
}

# Get all .fa (FASTA) files in the subdirectories of the base directory
$faFiles = Get-ChildItem -Path $baseDir -Recurse -Filter "*.fa"

# Initialize the counter for the number of .fa files found
$faFileCount = 0
Write-Output "----------------------------------- Step 1. Run antiSMASH for FASTA files -------------------------"
# Loop through each .fa file
foreach ($faFile in $faFiles) {
    # Increment the counter for each .fa file
    $faFileCount++

    # Get the directory of the current .fa file
    $fileDir = $faFile.DirectoryName

    # Get the folder name where the file is located and extract the part before the second underscore
    # This is done to reduce the name of the ouput folder
    $folderName = [System.IO.Path]::GetFileName($fileDir)
    $outputFolderName = $folderName.Split("_")[0] + "_" + $folderName.Split("_")[1]

    # Define the output subdirectory based on the extracted part of the folder name
    $outputSubDir = Join-Path -Path $outputDir -ChildPath $outputFolderName

    # Create the subdirectory for this file if it does not exist
    if (-Not (Test-Path -Path $outputSubDir)) {
        New-Item -ItemType Directory -Path $outputSubDir
    }


    # Output the file path information for testing
    Write-Output "------------------------------------------------------------------------------------------------"
    Write-Output " Processing FASTA file: $($faFile.FullName)"
    Write-Output " Output will be located in folder: $outputSubDir"
    Write-Output "------------------------------------------------------------------------------------------------"

    # Construct the Docker command as a single string
    $dockerCommand = "docker run --rm --volume ""${fileDir}:/input"" --volume ""${outputSubDir}:/output"" " +
                 "antismash/standalone $($faFile.Name) --output-dir /output " +
                 "--genefinding-tool prodigal-m "

    # Execute the Docker command
     Invoke-Expression $dockerCommand

    $optimise = "docker system prune -f "
     Invoke-Expression $optimise
}

# Output the total number of .fa files found
Write-Output "-------------------------------- antiSMASH done for all FASTA files ----------------------------------"
Write-Output "Total number of FASTA files found: $faFileCount"

# Python script to create plots and csv files with relevant information
# $targetDirectory = "C:\Users\a..."
$subdirectories = Get-ChildItem -Path $targetDirectory -Directory

# Loop through each subdirectory
foreach ($subDir in $subdirectories) {

    # Get all .gbk files in the current subdirectory
    $nodeGbkFiles = Get-ChildItem -Path $subDir.FullName -Recurse -Filter "*.gbk"

    # Check if any files were found in the current subdirectory
    if ($nodeGbkFiles.Count -eq 0) {
        Write-Output "No .gbk files starting with 'NODE' were found in $($subDir.FullName)."
    } else {
        Write-Output "Found NODE files in directory: $($subDir.FullName)"

        # Create the 'Nodes' folder inside the subdirectory if it doesn't already exist
        $nodesFolder = Join-Path -Path $subDir.FullName -ChildPath "Nodes"
        if (-Not (Test-Path -Path $nodesFolder)) {
            New-Item -ItemType Directory -Path $nodesFolder
            Write-Output "Created folder: $nodesFolder"
        }

        # Move each NODE .gbk file into the 'Nodes' folder
        foreach ($file in $nodeGbkFiles) {
            $destinationPath = Join-Path -Path $nodesFolder -ChildPath $file.Name
            Move-Item -Path $file.FullName -Destination $destinationPath -Force
        }


        # After processing, run Docker command inside each Nodes folder
        Write-Output "Running Docker command in the Nodes folder: $nodesFolder"

        # Define variables for the input and output directories for Docker
        $inputDir = $nodesFolder  # Use the Nodes folder as input
        $outputDir = Join-Path -Path $nodesFolder -ChildPath "Trees"  # Output will be in a folder inside Nodes

        # Create the output directory if it doesn't exist
        if (-Not (Test-Path -Path $outputDir)) {
            New-Item -ItemType Directory -Path $outputDir
        }

        # Construct the Docker command
        $dockerCommand = "docker run --rm " +
                         "--volume `"" + $inputDir + ":/home/input`" " +
                         "--volume `"" + $outputDir + ":/home/output`" " +
                         "currocam/big-scape:v1.1.5 " +
                         "--input /home/input --output /home/output --cutoffs 0.3 --include_singletons"

        # Execute the Docker command
        Invoke-Expression $dockerCommand

    }
}
Write-Output "-------------------------------- BiG-SCAPE done ------------------------------------"
