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

