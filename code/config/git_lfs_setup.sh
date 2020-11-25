#!/bin/bash

# This script sets up Git LFS to handle large files.

# Following the instructions in the docs here: https://git-lfs.github.com/

# Install Git LFS for the
git lfs install

# Track specific file types that are usually large (can add more)
git lfs track "*.rds" "*.ipynb"

# Make sure .gitattributes is tracked
git add .gitattributes
