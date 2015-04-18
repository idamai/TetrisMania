#!/bin/bash

MAXNODES=$1

###############
## BOOTSTRAP ##
###############
pm2 delete all
mvn install

## Start nodes ########################
for i in $(eval echo {1..$MAXNODES})
do
  NAME="NODE${i}"
  echo "Starting ${NAME}.."
  pm2 start src/main/dev/startNode.sh --name "$NAME" -f -- $NAME
done
## /Start nodes ########################

