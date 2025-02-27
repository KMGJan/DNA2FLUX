#!/bin/bash

# Check if a directory argument is provided
if [ -z "$1" ]; then
  echo "No directory specified. Using current directory."
  output_dir="."
else
  output_dir="$1"
fi

url="https://shark.smhi.se/api/sample/dataset"

# Perform the POST request to get the dataset keys
dataset_keys=$(curl -X 'POST' \
  "$url" \
  -H 'accept: application/json' \
  -H 'Content-Type: application/json' \
  -d '{
  "query": {
    "bounds": [],
    "fromYear": 2007,
    "toYear": 2030,
    "months": [],
    "dataTypes": [
      "Zooplankton",
      "Phytoplankton"
    ],
    "parameters": [],
    "checkStatus": "",
    "qualityFlags": [],
    "deliverers": [],
    "orderers": [],
    "projects": [],
    "datasets": [],
    "minSamplingDepth": "",
    "maxSamplingDepth": "",
    "redListedCategory": [],
    "taxonName": [],
    "stationName": [
      "BY31 LANDSORTSDJ"
    ],
    "vattenDistrikt": [],
    "seaBasins": [],
    "counties": [],
    "municipalities": [],
    "waterCategories": [],
    "typOmraden": [],
    "helcomOspar": [],
    "seaAreas": []
  }
}' | jq -r '.datasets[].key')

# Loop through the dataset keys and download each dataset
for dataset in $dataset_keys; do
  out="${dataset%.zip}"

  # Skip files containing "Algtox" in the name
  if [[ "$out" == *"Algtox"* ]]; then
    echo "Skipping file $out (contains 'Algtox')"
    continue
  fi

  file_path="$output_dir/$out.csv"

  # Extract components from the dataset key
  IFS='_' read -r prefix datatype year institute version date <<< "$out"

  # Construct the base name for non-PhysicalChemical format: SHARK_DataType_Year_INSTITUTE
  base_name="${prefix}_${datatype}_${year}_${institute}"

  # Check if an older version already exists
  existing_files=($(ls "$output_dir"/"$base_name"_*.csv 2>/dev/null))

  if [[ ${#existing_files[@]} -gt 0 ]]; then
    echo "Existing files found for $base_name:"

    # Find the latest version based on the filename
    latest_existing_file=$(ls "$output_dir"/"$base_name"_*.csv | sort -r | head -n1)

    # Extract the date part of the latest existing file
    latest_existing_date=$(echo "$latest_existing_file" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}')

    # Extract the date from the new file
    new_file_date=$(echo "$out" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}')

    echo "Existing: $latest_existing_date"
    echo "Newest: $new_file_date"

    # Compare dates and replace old files if a newer version is found
    if [[ "$new_file_date" > "$latest_existing_date" ]]; then
      echo "Newer version found. Replacing old file..."
      rm -f "${existing_files[@]}"
    else
      echo "No newer version available. Skipping download..."
      continue
    fi
  fi

  # Download the new file
  download_url="$url/$dataset/table"
  echo "Downloading: $download_url"
  curl -X 'GET' "$download_url" -H 'accept: application/json' | jq -r '(.headers | @csv), (.rows[] | @csv)' > "$file_path"
  echo "Saved as $file_path"

done
