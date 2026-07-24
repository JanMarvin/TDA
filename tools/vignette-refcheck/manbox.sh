#!/bin/sh
# usage: manbox.sh SECTION [lines]
s=$1; n=${2:-160}
start=$(grep -n "^$s " /home/claude/man/manual.txt | tail -1 | cut -d: -f1)
sed -n "$start,$((start+n))p" /home/claude/man/manual.txt | grep -v '^\s*$' | grep -v '\.tex\|^\f' | cut -c1-100
