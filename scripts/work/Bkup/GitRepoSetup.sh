#!/bin/bash

# 24/02/21 - Stephen JD Kay - University of Regina
# This script sets up the relevant git repos for use in analysing the pion-lt data
# It assumes you have hallc_replay_lt, UTIL_PION and UTIL_BATCH forked to your personal github account
# Input for the script is simply your github username

### 17/03/21 - SK - NOTE - You may not have the offline_pionlt in your fork, as such, you may need to get this first before running the script
### Some basic instructions on how to get a branch you don't have below
### git clone YOUR_FORK_OF_REPO ### Use the usual here, just insert (in "") your repo on github
### cd YOUR_REPO
### git checkout -b offline_kaonlt
### git remote add upstream "https://github.com/JeffersonLab/REPO"
### git fetch upstream
### git merge upstream/offline_kaonlt
### git push -u origin offline_kaonlt
### You'll have to do this for each repo you want unfortunately, I would just add the new branch to all of them and then delete them after, I'll add this as an automated process to this script soon

### Note for myself on this - git ls-remote --heads https://github.com/REPO BRANCH | wc -l
### Returns 0 if branch doesn't exist, 1 if it does
### Added this in now and testing it

echo "Beginning setup of git repos in /group/c-kaonlt/USERS"
GITUSER=$1
if [[ -z "$1" ]]; then 
    echo "I need a GitHub username as an argument, please re-run with this as an argument"
    exit 2
fi
echo "GitHub username taken as ${GITUSER}"
# Set path depending upon hostname. Change or add more as needed  
if [[ "${HOSTNAME}" = *"farm"* ]]; then 
    GROUPPATH="/group/c-kaonlt/USERS"
    # If user has not yet made their own user directory in the USERS area, then make it
    if [ ! -d "${GROUPPATH}/${USER}" ]; then
	mkdir "${GROUPPATH}/${USER}"
    fi
    USERPATH="${GROUPPATH}/${USER}"
    cd "${USERPATH}"
else echo "Host not recognised, please add relevant pathing for hostname to the script and re-run"
fi
# Defines repos to clone, add more or change as needed, these should all be forked to your account
REPLAY_REPO="https://github.com/${GITUSER}/hallc_replay_lt"
UTIL_PION_REPO="https://github.com/${GITUSER}/UTIL_PION"
UTIL_KAON_REPO="https://github.com/${GITUSER}/UTIL_KAONLT"
UTIL_BATCH_REPO="https://github.com/${GITUSER}/UTIL_BATCH"

# These are the active branches for each repo, change them as needed!
REPLAY_ACTIVE="offline_kaonlt"
PION_ACTIVE="offline_kaonlt"
KAON_ACTIVE="offline_kaonlt"
BATCH_ACTIVE="offline_kaonlt"
# Info dump to screen for user to check and confirm
echo ""
echo "Attempting to clone and setup the following repositories and branches -"
echo "${REPLAY_REPO} - ${REPLAY_ACTIVE}"
echo "${UTIL_PION_REPO} - ${PION_ACTIVE}"
echo "${UTIL_KAON_REPO} - ${KAON_ACTIVE}"
echo "${UTIL_BATCH_REPO} - ${BATCH_ACTIVE}"
echo "Ensure these exists and that these are the active branches you want!"
echo""

# Confirmation prompt for user to check the above printout is correct before continuing
read -p "Proceed with setup? <Y/N> " prompt
if [[ $prompt == "y" || $prompt == "Y" || $prompt == "yes" || $prompt == "Yes" ]]; then
    echo "Cloning and setting up branches and upstream remotes as specified above"
else 
    echo "Please update script with desired changes to repos/branches and re-run if desired"
    exit 1
fi
# If the directory doesn't already exists, clone it
if [ ! -d "${USERPATH}/hallc_replay_lt" ]; then
    git clone "${REPLAY_REPO}"
else echo "hallc_replay_lt already exists, assuming this is correct and continuing setup"
fi
cd "${USERPATH}/hallc_replay_lt"
# Switch to active branch
# This if statement checks if your fork has the desired branch, if it does, it just checks it out and adds the remote
# If you do NOT have the branch, it gets it from upstream and adds it to your fork
if [ $(git ls-remote --heads ${REPLAY_REPO} ${REPLAY_ACTIVE} | wc -l) = 1 ]; then
    git checkout "${REPLAY_ACTIVE}"
    git remote add "upstream" "https://github.com/JeffersonLab/hallc_replay_lt"
elif [ $(git ls-remote --heads ${REPLAY_REPO} ${REPLAY_ACTIVE} | wc -l) = 0 ]; then
    echo "${REPLAY_ACTIVE} Branch does not exist in your fork of hallc_replay_lt, grabbing it"
    git checkout -b ${REPLAY_ACTIVE}
    git remote add "upstream" "https://github.com/JeffersonLab/hallc_replay_lt"
    git fetch upstream
    git merge upstream/${REPLAY_ACTIVE}
    git push -u origin ${REPLAY_ACTIVE}
fi

if [ ! -d "${USERPATH}/hallc_replay_lt/UTIL_PION" ]; then
    git clone "${UTIL_PION_REPO}"
else echo "UTIL_PION already exists, assuming this is correct and continuing setup"
fi
cd "${USERPATH}/hallc_replay_lt/UTIL_PION"
if [ $(git ls-remote --heads ${UTIL_PION_REPO} ${PION_ACTIVE} | wc -l) = 1 ]; then
    git checkout "${PION_ACTIVE}"
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_PION"
elif [ $(git ls-remote --heads ${UTIL_PION_REPO} ${PION_ACTIVE} | wc -l) = 0 ]; then
    echo "${REPLAY_ACTIVE} Branch does not exist in your fork of UTIL_PION, grabbing it"
    git checkout -b ${PION_ACTIVE}
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_PION"
    git fetch upstream
    git merge upstream/${PION_ACTIVE}
    git push -u origin ${PION_ACTIVE}
fi
cd "${USERPATH}/hallc_replay_lt"

if [ ! -d "${USERPATH}/hallc_replay_lt/UTIL_KAONLT" ]; then
    git clone "${UTIL_KAONLT_REPO}"
else echo "UTIL_KAONLT already exists, assuming this is correct and continuing setup"
fi
cd "${USERPATH}/hallc_replay_lt/UTIL_KAONLT"
if [ $(git ls-remote --heads ${UTIL_KAON_REPO} ${KAON_ACTIVE} | wc -l) = 1 ]; then
    git checkout "${KAON_ACTIVE}"
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_KAONLT"
elif [ $(git ls-remote --heads ${UTIL_KAON_REPO} ${KAON_ACTIVE} | wc -l) = 0 ]; then
    echo "${REPLAY_ACTIVE} Branch does not exist in your fork of UTIL_KAONLT, grabbing it"
    git checkout -b ${KAON_ACTIVE}
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_KAONLT"
    git fetch upstream
    git merge upstream/${KAON_ACTIVE}
    git push -u origin ${KAON_ACTIVE}
fi
cd "${USERPATH}/hallc_replay_lt"

if [ ! -d "${USERPATH}/hallc_replay_lt/UTIL_BATCH" ]; then
    git clone "${UTIL_BATCH_REPO}"
else echo "UTIL_BATCH already exists, assuming this is correct and continuing setup"
fi
cd "${USERPATH}/hallc_replay_lt/UTIL_BATCH"
if [ $(git ls-remote --heads ${UTIL_BATCH_REPO} ${BATCH_ACTIVE} | wc -l) = 1 ]; then
    git checkout "${BATCH_ACTIVE}"
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_BATCH"
elif [ $(git ls-remote --heads ${UTIL_BATCH_REPO} ${BATCH_ACTIVE} | wc -l) = 0 ]; then
    echo "${REPLAY_ACTIVE} Branch does not exist in your fork of UTIL_BATCH, grabbing it"
    git checkout -b ${BATCH_ACTIVE}
    git remote add  "upstream" "https://github.com/JeffersonLab/UTIL_BATCH"
    git fetch upstream
    git merge upstream/${BATCH_ACTIVE}
    git push -u origin ${BATCH_ACTIVE}
fi
cd "${USERPATH}/hallc_replay_lt"

echo "Switching c-pionlt pathing in scripts to c-kaonlt (if it exists)"

find "${USERPATH}/hallc_replay_lt/" -name "*.py" -print0 | xargs -0 sed -i "s/c-pionlt/c-kaonlt/"
find "${USERPATH}/hallc_replay_lt/" -name "*.sh" -print0 | xargs -0 sed -i "s/c-pionlt/c-kaonlt/" 
find "${USERPATH}/hallc_replay_lt/" -name "*.C" -print0 | xargs -0 sed -i "s/c-pionlt/c-kaonlt/" 

echo ""
echo "Setup finished, you should ensure that your personal fork is up to date"
echo ""
echo "To update your fork, follow instructions from https://hallcweb.jlab.org/wiki/index.php/ROOT_Analyzer/Git - Note you should use the relevant active branch and not develop as stated"
echo ""
echo "For resolving merge conflicts, you may find the following commands useful -"
echo "git checkout --theirs PATH_TO_FILE"
echo "git checkout --ours PATH_TO_FILE"
echo "These override the conflicting file with their file or your file respectively"

exit 0
