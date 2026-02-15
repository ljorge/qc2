#!/bin/bash

clear
make clean    && \
make          && \
make cppcheck && \
make coverage && \
make doc
