#!/bin/bash
#Original taken from: http://pebblesinthesand.wordpress.com/2008/05/22/a-srcipt-for-running-processes-in-parallel-in-bash/
NUM=0
QUEUE=""
MAX_NPROC=2 # default
REPLACE_CMD=0 # no replacement by default
USAGE="A simple wrapper for running processes in parallel. Takes as second argument the name of the file containing the set of commands to run
Usage: `basename $0` [-h] [-j nb_jobs] command arg_list
       -h	     Shows this help
       -j nb_jobs    Set number of simultanious jobs [2]
 Example:
 	`basename $0` -j 3 -f parallelcommands.tmp
"

function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}

function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
	if [ -d /proc/$PID  ] ; then
	    QUEUE="$QUEUE $PID"
	    NUM=$(($NUM+1))
	fi
    done
}

function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
	if [ ! -d /proc/$PID ] ; then
	    regeneratequeue # at least one PID has finished
	    break
	fi
    done
}

# parse command line
if [ $# -eq 0 ]; then #  must be at least one arg
    echo "$USAGE" >&2
    exit 1
fi

while getopts j:f:h OPT; do # "j:" waits for an argument "h" doesnt
    case $OPT in
	h)	 echo "$USAGE"
    	    exit 0 ;;
	j)   MAX_NPROC=$OPTARG ;;
	f)   COMMANDSFILE=$OPTARG ;;
	\?)  # getopts issues an error message
	    echo "$USAGE" >&2
	    exit 1 ;;
    esac
done

# Main program
echo Using $MAX_NPROC parallel threads
   
exec<$COMMANDSFILE

while read line
do
    #echo $line
    CMD=$line

    echo "Running $CMD" 
    
    #`$CMD` &
    eval $CMD &
    
    PID=$!
    queue $PID
    
    while [ $NUM -ge $MAX_NPROC ]; do
	checkqueue
	sleep 0.4
    done
done
wait # wait for all processes to finish before exit
