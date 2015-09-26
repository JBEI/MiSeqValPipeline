#!/bin/sh

curl -H "Content-Type: application/json" \
  -X POST \
  -d  '{email: "oge", password: "scot!syeP6"}' \
  https://registry.jbei.org/rest/accesstoken | jq '.sessionId'
