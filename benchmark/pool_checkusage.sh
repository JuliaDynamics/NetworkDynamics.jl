#!/bin/bash

# Variables
REMOTE_USER="wuerfelh"
REMOTE_HOST="$1.physik.hu-berlin.de"
THRESHOLD=10    # CPU usage threshold for reporting overall usage
DURATION=10     # Duration to monitor CPU usage in seconds
SINGLE_CORE_THRESHOLD=10 # CPU usage threshold for single process

# SSH into the remote machine and run the commands
ssh $REMOTE_USER@$REMOTE_HOST << EOF
  echo "Checking CPU usage on $REMOTE_HOST for $DURATION seconds..."

  # Capture top output for a specified duration
  top -b -d 2 -n $((DURATION / 2)) > top_output.txt

  echo "Checking for high overall CPU usage by user processes..."
  echo ""

  # Calculate average CPU usage for user processes
  AVG_CPU=\$(grep '^%Cpu' top_output.txt | awk '{ sum += \$2 } END { if (NR > 0) print sum / NR }')

  echo "Average CPU usage by user processes: \$AVG_CPU%"

  if (( \$(echo "\$AVG_CPU > $THRESHOLD" | bc -l) )); then
    echo "High overall CPU usage detected. Reporting processes:"
    ps -eo user,pid,ppid,%cpu,cmd --sort=-%cpu | awk -v threshold=$THRESHOLD '\$4 > threshold && \$1 != "'$REMOTE_USER'" { print }'
  else
    echo "No significant overall CPU usage detected."
  fi

  echo ""
  echo "Checking for any single process using $SINGLE_CORE_THRESHOLD % CPU..."
  echo ""

  # Check for any single process using more than SINGLE_CORE_THRESHOLD
  grep -E '^ *[0-9]+' top_output.txt | awk '{ if (\$9 >= $SINGLE_CORE_THRESHOLD) print \$1, \$2, \$9, \$12 }' | while read -r user pid cpu cmd
  do
    if [ "\$user" != "$REMOTE_USER" ]; then
      echo "Process \$pid (\$cmd) by user \$user is using \$cpu% CPU."
    fi
  done

  # Clean up
  rm top_output.txt
EOF
