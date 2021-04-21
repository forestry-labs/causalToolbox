#!/bin/bash
R CMD BATCH --no-save --args -$1 simulate.R  simulate$1.out
