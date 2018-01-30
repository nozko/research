#!/bin/sh

cat learn.sh | xargs -P7 -I{} -t bash -c '{}'
