#!/usr/bin/env bash
emcc smallpt.cpp --js-library library.js -s WASM=1 -O3 -o smallpt.html 