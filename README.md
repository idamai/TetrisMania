CS3243 Tetris Mania
=====================================

## Info

Uses maven for dependencies, junit testing, JAR executable deployment and creation.


## Running unit tests

    $ mvn test


## Creating JAR executable

    $ mvn install
    

## Steps to getting into development

1. Fork the project on Github.
2. Clone the fork
3. Link it back to main repository `upstream`

```
    $ git remote add upstream git@github.com:idamai/TetrisMania.git
```
    
    
4. To retrieve updates from main repository daily, use the following commands


```
    $ git fetch upstream              ## fetch changes from upstream
    $ git checkout master             ## checkout master branch
    $ git rebase upstream master      ## Syncs changes
```
   
    
5. Contribute to development on separate branch @ `origin/<branch name>` ; merge with origin when done.
6. Submit pull request.
7. Download Maven.



## Running the batch script

Batch script uses `pm2` and `nodejs`.  To install,

    $ npm install -g pm2

Run the script:

    ## 10 is the number of nodes to run
    $ ./dev 10

See the running processes:

    $ pm2 l
    
See the logs

    $ pm2 logs



## Distributed computing

In addition to normal GA, our algorithm features a distributed computing algorithm over multiple
computers, following a Master-Slave model.  This is done via sending UDP packets across the network.
The best weights are then gathered and run.

### Running slave nodes

All slaves must be run before the master is run.  This is to allow the notifications to be propagated to master.

    $ java -cp . -jar TetrisMania-1.0-SNAPSHOT-jar-with-dependencies.jar slave <HOST IP> <HOST PORT> <SLAVE PORT>
    $ ## Example
    $ java -cp . -jar TetrisMania-1.0-SNAPSHOT-jar-with-dependencies.jar slave localhost 9000 9001


### Running master

Master then runs.  It pends until all available slaves have submitted their results.

    $ java -jar TetrisMania-1.0-SNAPSHOT-jar-with-dependencies.jar master <INPUT PORT>
    $ ## Example
    $ java -jar TetrisMania-1.0-SNAPSHOT-jar-with-dependencies.jar master 9000


