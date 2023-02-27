#!/bin/bash

from_clarity=true

if [ "$from_clarity" = true ]; then
	rsync -v -r spirit:/home/bfildier/analyses/explore-multiscale-extremes/movies/ ../movies/
fi

exit 0
