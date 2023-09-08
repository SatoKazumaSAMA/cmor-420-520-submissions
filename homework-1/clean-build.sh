#!/bin/bash

if [ ! -d .build]

then
    mkdir .build

fi

pdflatex $1.tex

mv $1.aux $!1.log .build
