#!/bin/bash

NODE_NAME=$1
DIR_NAME=${PWD}
echo "Started $NODE_NAME."


##################################
## Insert program execution here

java -jar ${PWD}/target/TetrisMania-1.0-SNAPSHOT-jar-with-dependencies.jar

## End program execution here
##################################



## Stop node once done
pm2 stop "$NODE_NAME"
