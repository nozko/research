#!/bin/sh

cat learn.sh | xargs -P4 -I{} -t bash -c '{}'
