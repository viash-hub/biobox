#!/bin/bash

arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s\(.*\)/arriba: "\1"/'