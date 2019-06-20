#!/bin/bash

echo $1
hg commit -m "$1"
hg push ssh://hg@bitbucket.org/LMalavolta/harps_skycorr
