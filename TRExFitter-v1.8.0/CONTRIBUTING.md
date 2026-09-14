# Contributing to TRExFitter

## Overview
1. [Just show me the commands](#just-show-me-the-commands)
2. [How to contribute](#how-to-contribute)
3. [Code style](#code-style)
    - [General](#general)
    - [Style](#style)
    - [Safety](#safety)
4. [Git how to](#git-how-to)
    - [Setup git](#setup-git)
    - [Fork TRExFitter repository on gitlab](#fork-trexfitter-repository-on-gitlab)
    - [Setup a clean development area](#setup-a-clean-development-area)
    - [Commit changes to your local repository and push them to gitlab](#commit-changes-to-your-local-repository-and-push-them-to-gitlab)
    - [Create merge request](#create-merge-request)
    - [Other useful commands](#other-useful-commands)
    - [Helpful links](#helpful-links)
5. [Keeping track of the changes](#keeping-track-of-the-changes)
6. [Continuous integration](#continuous-integration)
    - [Static code analysis](#static-code-analysis)
    - [Tests of the framework](#tests-of-the-framework)
    - [Creation of docker images](#creation-of-docker-images)

## Just show me the commands
For a quick svn-git translation table, see the [git-svn crash course](https://git.wiki.kernel.org/index.php/GitSvnCrashCourse).

## How to contribute
This doc is mainly for developers, hopefully serving as a quick reference for using git.
GitLab repos can be used in numerous ways, but the recommended way is to set up your GitLab account with ssh keys and use the ssh protocol.

## Code style
Please try to follow the rules below for your merge request.

### General
* Make the code modular - try to split the code into smaller functions can be reused multiple times (smaller function means < 100 lines)
* Avoid using pointers (allocating on heap) - pointers should be used in cases where they can point to null - e.g. when reading something that may not exist
* Do not ignore compilation warnings - fix them
* Use const as much as possible - in function arguments, when declaring variables, etc
* Functions should be marked as const when appropriate
* When function doesn’t change its arguments they should be passed as const reference
* Never use `using` in header space
* Try not to use `using namespace std` even in source files - this can cause problems
* Each header file needs to have an include guard
* Keep lines of reasonable length
* Use forward class declarations for header files whenever possible
* Avoid using raw pointers, use smart pointers
* Do not call destructors directly, e.g. `h->~TH1D()` (this _does not_ release the memory!), use standard `delete h`

### Style
* Use 4 spaces for indentation - do not use tabs
* Classes start with capital letters
* Member variables start with a lowercase “f” (this is not very common but this was used in TRExFitter and we should stick to it) then follow by an uppercase letter
* Use “nullptr” for null pointer values, do not use `0` , `0x0` or `NULL`
* Use “{}” even for one-line statement after if/else
* Remove code if it is not needed, do not just comment it out
* Avoid using `or` in `if` statements, use `||` and `&&` instead of `and`
* `else` should be on next line after `if`
* Use doxygen documentation style for comments in headers
* Try to keep the style of variable names used in the code - generally using CamelCase but with variables ending with `up`/`down` it is ok to use `_`, example: `TotalSyst_up`
* Do not leave trailing-whitespaces in the code

### Safety
* Do not let the code crash, always protect against possible failures
* Limit variable scope - define variable only when it is needed

## Git how to
### Setup git
On lxplus, git is already installed.
Just check what version you are using:
```bash
git --version
```

For first time use, you must setup your configuration, [cern gitlab "getting started"](https://cern.service-now.com/service-portal/article.do?n=KB0003137).
In addition to username and email might be interested in setting your preferred editor and colored diff.
However, don't bother with the password caching.
Instead, [setup your ssh keys](https://gitlab.cern.ch/help/ssh/README.md) or use kerberos.
```bash
git config --global user.name "First Last"
git config --global user.email "first.last@cern.ch"
git config --global color.ui auto
git config --global core.editor "emacs -nw"
```
You might prefer `vim` as an editor, in which case you update the above accordingly.

You might have to require us to add you as a developer before you can push changes to the gitlab repository.
In the meantime, you can still commit them to your local repository.

### Fork TRExFitter repository on gitlab
Click on the `Fork` button on the main gitlab page.

### Setup a clean development area
```bash
git clone ssh://git@gitlab.cern.ch:7999/<yourname>/TRExFitter.git
```

### Commit changes to your local repository and push them to gitlab
```bash
cd TRExFitter
git status # start from a clean state, up-to-date with the origin
git remote add upstream ssh://git@gitlab.cern.ch:7999/TRExStats/TRExFitter.git # sets upstream to the main repository, you only need to do this once
git fech upstream # downloads the changes from upstream repository
git rebase upstream/master # rebase on top of the upstream/master branch
git checkout -b some-project  # put your work on a branch
git diff <file> # look at the diff between current HEAD and local version to see what you will add
git add <list of files> # add files you want to change
git commit -m "your commit message" some.txt file.txt  # commit files and prompt editor for commit message
git push -u origin some-project  # push your changes to gitlab, parameter -u sets upstream
git checkout master  # to go back to the origin, to be able to synchronise with it again
```

### Create merge request
When your branch is pushed to your fork (origin).
You should be able to make a merge request into the upstream.
Just click on `Create merge request` button and add description.

### Other useful commands:
```bash
git pull # is combination of git fetch and git merge, example: git pull upstream master
git checkout -- <files> # remove local changes on files in <files>
git merge <branch> # merges <branch> into your current branch
git branch -a # prints all local and remote branches
git branch -D <branch> # delete local branch <branch>
```

Notes:
- Feature development should always occur on a dedicated branch, rather than on the master branch.
- Try to pick a branch name that is short and describes what you are trying to do.
- When you have made enough progress, push to gitlab, then issue a merge request (green button), and start discussing the proposed changes.
- Recommendations about commit messages: Follow Linus' guidelines: one short descriptive line, followed by a longer description of the motivation and implementation, if necessary. Limit line length (with emacs `M-q` or `M-x fill-paragraph`).
- Setting up your git preferences (editor, aliases, etc.): see the [Git Configuration](http://git-scm.com/book/en/Customizing-Git-Git-Configuration) docs.
- Report issues and bugs: open a new issue, on gitlab or JIRA

### Helpful links
- [Reference of all git commands](http://git-scm.com/docs)
- When you have time, I recommend reading through [the git book](http://git-scm.com/book).
- [Here's a simple guide for using git](http://rogerdudler.github.io/git-guide/), found by Suneet.
- For more details on the workflow, see [a simple git branching model](https://gist.github.com/jbenet/ee6c9ac48068889b0912) or [a detailed description of the feature-branch](https://www.atlassian.com/git/workflows#!workflow-feature-branch).

## Keeping track of the changes
This project is linked to the [TRExFitter JIRA project](https://its.cern.ch/jira/projects/TTHFITTER).
You can interact with the JIRA project by quoting the ticket ID in your commit messages or merge requests message.
```bash
git commit -m "Working on TTHFITTER-102 - Fixing some pointer initialisation"
```

This will automatically create a comment on the corresponding JIRA ticket.
Some keywords are also usable to close the JIRA ticket when needed.
Those keywords are:
```bash
Closes TTHFITTER-102
Resolves TTHFITTER-102
Fixes TTHFITTER-102
```

Such messages can be put in the commit messages or in the merge request message.

This action is only performed when the users branch has been merged to the master.
**It is strongly recommended to use this feature as much as possible to make easier the tracking of the various changes.**

## Continuous integration
TRExFitter automatically runs a continuous integration pipeline when you commit changes to the repository.
The steps are defined in `.gitlab-ci.yml`.
There are currently three different types of operations:
- static code analysis
- tests of the framework
- creation of docker images

### Static code analysis
Two types of static code analysis checks are performed, using [CPD](https://pmd.github.io/latest/pmd_userdocs_cpd.html) and [cppcheck](http://cppcheck.sourceforge.net/).
The number of issues flagged is compared to a reference, and these CI stages fail if the number changes.
New merge requests should avoid introducing new issues.
If they fix existing issues, the reference numbers can be updated in the scripts found in `test/scripts/`.

### Tests of the framework
Tests of the framework are automatically performed to spot possibly bugs introduced by accident.
This includes compilation and multiple example workflows with different fit models and options.
The relevant files are found in `test/`.
Occasionally it is possible that an update causes checks to fail, even if there is nothing wrong with the update.
An example for this is switching the ROOT version, which might cause some differences due to finite floating point accuracy in calculations.
The script `test/scripts/dump_log_files.sh` can be used to update the reference output expected in the CI.
Any change in reference files should be carefully justified.

### Creation of docker images
A docker image tagged `latest` is created to correspond to the current state of the master branch.
Images are also created for tagged version of TRExFitter.
You can find all images in the [container registry](https://gitlab.cern.ch/TRExStats/TRExFitter/container_registry).