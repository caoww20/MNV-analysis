#!/bin/bash

# Check number of arguments
if (( $# < 2 )); then
  echo "Please provide a command and a second string parameter"
  exit 1
fi
# Remove data
# rm -r all_res2 mem.txt MNVAnno_time.txt time.txt

# Record start time
start_time=$(date +%s%N)

# Read arguments: first argument is the command, second argument is the ps search string
command="$1"
string="$2"

# Record start time
start_time=$(date +%s%N)

# Run the command
# (time $command ) 2>>MNVAnno_time.txt &
$command &

# Record the process ID
pid=$!

# Initialize maximum memory usage
mem_max=0

# Check every second if the process has ended
while true; do
  if ! kill -0 $pid 2> /dev/null; then
    break
  fi
  # Use ps aux to get results for background processes containing the string (e.g., "VCFCheck" or "MNVIdentify"), extract the 6th column, sort, and get the last value
  # mem=$(ps aux | grep jinww|grep -E "bcftools" | grep -v grep | awk '{print $6}' | sort -n | tail -n 1)
  mem=$(ps aux | grep jinww|grep -E "$string" | grep -v grep | awk '{print $6}' | sort -n | tail -n 1)
  # Compare with mem_max and replace if larger
  if (( mem > mem_max )); then
    mem_max=$mem
  fi
done

# Record maximum memory usage to file
echo "Max memory of process $pid: $mem_max KB" >> mem.txt

# Record end time
end_time=$(date +%s%N)

# Calculate command execution time
elapsed_time=$(( (end_time - start_time) / 1000000 ))

# Record command execution time
echo "Elapsed time of process $pid: $elapsed_time ms" >> time.txt