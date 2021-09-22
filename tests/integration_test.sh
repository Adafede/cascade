#!/usr/bin/env bash

if [ ! -d config ]; then
  echo "Sorry, you need to run that from where your config is."
  exit 1
fi

cp -R config/default config/params &&
echo "Real tests will be implemented"

