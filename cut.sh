#!/bin/bash
cat all_counts.txt | cut -f1,7- | sed 1d > counts_matirx.txt
